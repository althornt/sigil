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

##################
# Functions
##################

importMetaKallisto <- function(row){
  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general, LM22)

  # Add metadata to column
  df_metadata$data_source <- row[1] # add name of data source
  df_metadata$type <- row[4] # add rna-seq type (paired vs single)

  # Get paths to kallisto files
  res_path <- file.path(row[2], "kallisto_out",df_metadata$Run, "abundance.h5")

  return(list(
      "ls_kallsto_paths"=res_path,
      "metadata"=df_metadata,
      "ls_samples_run"=df_metadata$Run))
}

make_umap <- function(num_neighbor,meta_col) {

  set.seed(123)

  # Make color palette
  n <- length(unique(df_merged_metadata_lm22[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, df_merged_metadata_lm22)

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
  sampleTable <- df_merged_metadata_lm22 %>%
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


###################
# MAIN
###################
# Arguments
option_list <- list(

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "path to write outputs"),

  optparse::make_option(
    c("-m", "--manifest"),
    type = "character",
    default = NULL,
    help = "path to tsv manifest file with 3 columns: data_source, res_path, metadata_path, type")

  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/UMAPs/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs/"), recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/deseq2_outputs/"))){
  dir.create(file.path(opt$out_dir,"/deseq2_outputs/"), recursive = TRUE, showWarnings = TRUE)}

# Open files
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(head(df_manifest))

# Import and combine source metadata files
ls_kallisto_meta = apply(df_manifest, 1, importMetaKallisto)

# Split into kallisto and metadata files for each data set
ls_kallisto <- ls_mesa <- ls_sample_names <- c()
for (item in ls_kallisto_meta) {
     ls_kallisto <- append(ls_kallisto, item[1])
     ls_mesa <- append(ls_mesa, item[2])
     ls_sample_names <- append(ls_sample_names, item[3])
   }

# Merge metadata from each data source by rows
df_merged_metadata <- do.call("rbind", ls_mesa)
rownames(df_merged_metadata) <- c()

# List of samples with LM22 labels
ls_smpls_lm22 <- df_merged_metadata %>%
   dplyr::filter(LM22 != "") %>%
   dplyr::pull(Run)

# Remove samples without LM22 labels
df_merged_metadata_lm22 <- df_merged_metadata %>%
   dplyr::filter(LM22 != "")

# Make list of kallisto paths from each data source
ls_kallisto_paths <- unlist(ls_kallisto,use.names = FALSE,recursive = FALSE)
names(ls_kallisto_paths) <- unlist(ls_sample_names)

# Remove kallisto paths without LM22 labels
ls_kallisto_paths_lm22 <- ls_kallisto_paths[ls_smpls_lm22]

print(length(ls_kallisto_paths_lm22))

# Import counts with tximport
if(all(file.exists(ls_kallisto_paths_lm22)) != TRUE)
  stop("Error: missing kallisto h5 files.")

# Transcripts to gene, used in tximport
edb <- EnsDb.Hsapiens.v86
tx2gene = transcripts(edb, columns=c("tx_id", "gene_name"),return.type="DataFrame")

txi.kallisto <- tximport(ls_kallisto_paths_lm22, type = "kallisto",tx2gene = tx2gene,
                         txOut = FALSE,ignoreTxVersion=TRUE,
                         countsFromAbundance = "scaledTPM")
dat = txi.kallisto$counts

# Write combined TPM file and log2 TPM
write.csv(dat,
         file.path(opt$out_dir,"combined_kallisto_tpm.csv"),
         row.names = TRUE)
log2trans_dat <- as.data.frame(log2(dat +1))
write.csv(log2trans_dat,
         file.path(file.path(opt$out_dir,"combined_kallisto_log2tpm.csv")),
         row.names = TRUE)


##################
# UMAPs
##################

# # Drop genes with low variance.
# getVar <- apply(log2trans_dat_lm22[, -1], 1, var)
# param <- median(getVar)
# log2trans_dat_filt <- log2trans_dat_lm22[getVar > param & !is.na(getVar), ]
#
# # Transpose and format
# log2trans_dat_filt_t <- as.data.frame(t(log2trans_dat_filt)) %>%
#   dplyr::filter(!row.names(.) %in% c("gene"))
# rownames(log2trans_dat_filt_t) <- colnames(log2trans_dat_filt)
#
# # PCA.
# prcomp.out = as.data.frame(prcomp(log2trans_dat_filt_t, scale = F)$x)
# prcomp.out$Run = rownames(log2trans_dat_filt_t)
# prcomp.out.merge = merge(prcomp.out, y = log2trans_dat_lm22)
#
# # Making variations of UMAPs with different numbers of neighbors
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source")


##################
# DESEQ2 LM22
##################

# Run deseq2 on each LM22 cell type vs all others
if("LM22" %in% colnames(df_merged_metadata_lm22)){
  print("Running...")

  # Format DF to label samples(columns) with general and more specific labels
  df_sample_annotations <- df_merged_metadata_lm22 %>%
    dplyr::select(Run,sigil_general,data_source,LM22,type) %>%
    tibble::column_to_rownames("Run")

  print(head(df_sample_annotations))

  # Run DE 1 cell type vs all others
  DE_results <- sapply(
    unique(df_merged_metadata_lm22[["LM22"]])[2],
    runDE_1_vs_all,
    meta_col_to_use="LM22")

  # print("Run DE done____________________ ")
  #
  # # Run heatmap function on all cell_types and unnest the list of DEG
  # all_DEG_res <- unlist(lapply(
  #   unique(df_merged_metadata_lm22[["LM22"]]), list2heatmap,
  #   meta_col_to_use="LM22",results=DE_results),
  #   recursive = TRUE, use.names = FALSE)
  #
  # # Combine data from individual comparisons
  # df_all_DEG <- log2trans_dat %>%
  #   tibble::rownames_to_column("gene") %>%
  #   dplyr::filter(gene %in% all_DEG_res) %>%
  #   tibble::column_to_rownames("gene")
  # write.csv(df_all_DEGdf_all_DEG,
  #           paste0(opt$out_dir,"/deseq2_outputs/LM22_all_DEG.csv"))
  #
  # # Make heatmap with all DEG from sigil_general
  # pheatmap_combined_all_deg <- pheatmap(
  #   df_all_DEG,
  #   main = paste0("All DEG "),
  #   scale = "row",
  #   show_rownames = F,
  #   show_colnames = F,
  #   annotation_col = df_sample_annotations)
  #
  # save_pheatmap_pdf(pheatmap_combined_all_deg,
  #                   paste0(opt$out_dir,"/deseq2_outputs/LM22_all_DEG_heatmap.pdf"))
}
