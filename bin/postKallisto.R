#!/usr/bin/env Rscript

library(optparse)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(uwot)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(ensembldb)
library(DESeq2)
library(pheatmap)

##############
# Functions
##############

# UMAP function
make_umap <- function(num_neighbor,meta_col) {

  set.seed(123)

  # Make color palette
  n <- length(unique(metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, metadata)

  # Plot UMAP
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP Kallisto: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste("UMAPs/kallisto_PCA_UMAP",meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)

}

# Run Deseq2 for one cell type vs all other samples
runDE_1_vs_all <- function(meta_col_to_use, cell_type_val) {

  print(cell_type_val)

  # Make sample table to compare given cell type vs all other cell types
  sampleTable <- metadata %>%
    dplyr::select(Run, meta_col_to_use) %>%
    dplyr::mutate(condition = ifelse(get(meta_col_to_use) == cell_type_val, "main", "other"))

  # Run DESEQ2
  dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
  colData(dds)$condition<-factor(colData(dds)$condition, levels=c("other", "main"))
  dds_<-DESeq(dds)
  res<-results(dds_)
  res<-res[order(res$padj),]

  #Write output to file
  write.csv(as.data.frame(res),file.path(opt$out_dir,"deseq2_outputs",paste0(cell_type_val,".csv")))
  print("wrote DE output")

  df_res <- as.data.frame(res) %>%
    tibble::rownames_to_column(var = "gene")

  # Get top DEG with positive and negative log2FC
  DEG_up <- df_res %>%
    dplyr::filter(padj < .001) %>%
    dplyr::filter(log2FoldChange > 2) %>%
    dplyr::pull("gene")

  DEG_down <- df_res %>%
    dplyr::filter(padj < .001) %>%
    dplyr::filter(log2FoldChange < -2 ) %>%
    dplyr::pull("gene")

  print("DE func done")
  return(list("DEG_down"= DEG_down, "DEG_up"= DEG_up))
}

# Function to save pheatmaps to a file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Function to make a heatmap from a DEG list from 1 cell type
list2heatmap <- function(cell_type, meta_col_to_use, results){

  # UP DEG
  # up_cell_type <- log2trans_dat %>%
  #   dplyr::filter(row.names(log2trans_dat) %in% unlist(results["DEG_up",cell_type]))

  up_cell_type <- log2trans_dat %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(gene %in% unlist(results["DEG_up",cell_type])) %>%
    tibble::column_to_rownames("gene")

  print("in list2heatmap")
  print(head(up_cell_type))
  print(dim(up_cell_type))

  # # DOWN DEG
  # down_cell_type <- log2trans_dat %>%
  #   dplyr::filter(row.names(log2trans_dat) %in% unlist(results["DEG_down",cell_type]))

  # UP and DOWN DEG
  # up_down_cell_type <- log2trans_dat %>%
  #   dplyr::filter(row.names(log2trans_dat) %in% unlist(results["DEG_down",cell_type])  |
  #                   row.names(log2trans_dat) %in% unlist(results["DEG_up",cell_type]))

  # If at least 2 DEG in the condition , make heatmap
  if (nrow(up_cell_type) >= 2){

    up_deg_heatmap <- pheatmap(
      main = paste0(" UP DEG in ", cell_type),
      up_cell_type,
      scale = "row",
      show_rownames=F,
      show_colnames=F,
      annotation_col=df_sample_annotations)

    save_pheatmap_pdf(
      up_deg_heatmap,
      paste0(opt$out_dir,"/deseq2_outputs/",meta_col_to_use,"_", cell_type,"_UP_DEG_heatmap.pdf"))
  }

  return(list("DEG"=row.names(up_cell_type)))

}

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--kallisto_dir"),
    type = "character",
    default = NULL,
    help = "path to directory with kallisto outputs"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "path to write outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "path to metadata file"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directory
if (!dir.exists(opt$out_dir)) {dir.create(opt$out_dir)}

# Open files
metadata = read.csv(file = opt$metadata)
files <- file.path(opt$kallisto_dir, metadata$Run, "abundance.h5")

# Transcripts to gene, used in tximport
edb <- EnsDb.Hsapiens.v86
tx2gene = transcripts(edb, columns=c("tx_id", "gene_name"),return.type="DataFrame")

# Check all kallisto files exist
if(all(file.exists(files)) != TRUE)
  stop("Error: Missing kallisto abundance.h5 files.")

names(files)<- metadata$Run
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene,
                         txOut = FALSE,ignoreTxVersion=TRUE,
                         countsFromAbundance = "scaledTPM")

# Write combined TPM file and log2 TPM
dat = txi.kallisto$counts
write.csv(dat,
          file.path(opt$out_dir,"combined_kallisto_tpm.csv"),
          row.names = TRUE)
log2trans_dat <- as.data.frame(log2(dat +1))
write.csv(log2trans_dat,
          file.path(file.path(opt$out_dir,"combined_kallisto_log2tpm.csv")),
          row.names = TRUE)

#################
# UMAP
#################
# Make output directory
if (!dir.exists(paste0(opt$out_dir,"/UMAPs/"))){
  dir.create(paste0(opt$out_dir,"/UMAPs/"))
}
# Drop genes with low variance.
getVar <- apply(log2trans_dat[, -1], 1, var)
param <- median(getVar)
log2trans_dat_filt <- log2trans_dat[getVar > param & !is.na(getVar), ]

# PCA.
prcomp.out = as.data.frame(prcomp(as.data.frame(t(log2trans_dat_filt)), scale = F)$x)
prcomp.out$Run = rownames(prcomp.out)
prcomp.out.merge = merge(prcomp.out, y = metadata)

# Making variations of UMAPs with different numbers of neighbors
# only make UMAPs if enough samples
if(("sigil_general" %in% colnames(metadata)) & (ncol(dat)> 4)){
  lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general")
}

# Color by cell type
if(("sigil_cell_type_treatment" %in% colnames(metadata)) & (ncol(dat)> 4))
{
  lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_cell_type_treatment");
}

########################
# DESEQ2 sigil_general
########################


# Make output directory
if (!dir.exists(paste0(opt$out_dir,"/deseq2_outputs/"))){
  dir.create(paste0(opt$out_dir,"/deseq2_outputs/"))
}

# Run deseq2 on each cell type vs all others
if("sigil_general" %in% colnames(metadata)){

  # DF to label samples(columns) with general and more specific labels if they exist
  if("sigil_cell_type_treatment" %in% colnames(metadata)){
    df_sample_annotations <- metadata %>%
      dplyr::select(Run,sigil_general, sigil_cell_type_treatment) %>%
      tibble::column_to_rownames("Run")
  }
  else {
    df_sample_annotations <- metadata %>%
      dplyr::select(Run,sigil_general) %>%
      tibble::column_to_rownames("Run")
  }

  # Run DE 1 cell type vs all others
  sigil_general_results <- sapply(
    unique(metadata[["sigil_general"]]),
    runDE_1_vs_all,
    meta_col_to_use="sigil_general")

  print(sigil_general_results)

  # Run heatmap function on all cell_types and unnest the list of DEG
  all_DEG_res_sigil_general <- unlist(lapply(
    unique(metadata[["sigil_general"]]), list2heatmap,
    meta_col_to_use="sigil_general",results=sigil_general_results),
    recursive = TRUE, use.names = FALSE)

  # Combine data from individual comparisons
  all_DEG_data_sigil_general <- log2trans_dat %>%
    dplyr::filter(row.names(log2trans_dat) %in% all_DEG_res_sigil_general)
  write.csv(all_DEG_data_sigil_general,
            paste0(opt$out_dir,"/deseq2_outputs/sigil_general_all_DEG.csv"))

  # Make heatmap with all DEG from sigil_general
  pheatmap_combined_all_deg_sigil_general <- pheatmap(
    all_DEG_data_sigil_general,
    main = paste0("All DEG"),
    scale = "row",
    show_rownames = F,
    show_colnames = F,
    annotation_col = df_sample_annotations)

  save_pheatmap_pdf(pheatmap_combined_all_deg_sigil_general,
                    paste0(opt$out_dir,"/deseq2_outputs/sigil_general_all_DEG_heatmap.pdf"))
}
######################################
# DESEQ2 sigil_cell_type_treatment
######################################

# Run deseq2 on each cell type (treatment specific) vs all others
if("sigil_cell_type_treatment" %in% colnames(metadata)){

  # DF to label samples(columns) with general and more specific labels
  df_sample_annotations <- metadata %>%
    dplyr::select(Run,sigil_general, sigil_cell_type_treatment) %>%
    tibble::column_to_rownames("Run")

  # Run DE 1 cell type vs all others
  sigil_cell_type_treatment_results <- sapply(
    unique(metadata[["sigil_cell_type_treatment"]]),
    runDE_1_vs_all,
    meta_col_to_use="sigil_cell_type_treatment")

  # Run heatmap function on all cell_types and unnest the list of DEG
  all_DEG_res_cell_type_treatment <- unlist(lapply(
    unique(metadata[["sigil_cell_type_treatment"]]), list2heatmap,
    meta_col_to_use="sigil_cell_type_treatment",results=sigil_cell_type_treatment_results),
    recursive = TRUE, use.names = FALSE)

  # Combine data from individual comparisons
  all_DEG_cell_type_treatment <- log2trans_dat %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(gene %in% all_DEG_res_cell_type_treatment) %>%
    tibble::column_to_rownames("gene")
  write.csv(all_DEG_cell_type_treatment,
            paste0(opt$out_dir,"/deseq2_outputs/sigil_cell_type_treatment_all_DEG.csv"))

  # Make heatmap with all DEG from sigil_general
  pheatmap_combined_all_deg_sigil_cell_type_treatment <- pheatmap(
    all_DEG_cell_type_treatment,
    main = paste0("All DEG using treatment"),
    scale = "row",
    show_rownames = F,
    show_colnames = F,
    annotation_col = df_sample_annotations)

  save_pheatmap_pdf(pheatmap_combined_all_deg_sigil_general,
                    paste0(opt$out_dir,"/deseq2_outputs/sigil_cell_type_treatment_all_DEG_heatmap.pdf"))
}
