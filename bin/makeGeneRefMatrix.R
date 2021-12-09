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
library(limma)

##################
# Functions
##################
read_in_deseq2 <- function(file_path, topN) {
  df_res <- read.csv(file = file_path, header = TRUE)

  # Get cell type name from file path
  str_cell_type <- tools::file_path_sans_ext(basename(file_path))

  # Make volcano plot
  volcano <- ggplot(data=df_res, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point() + theme_minimal() +
    geom_hline(yintercept=-log10(0.05), col="red") +
    labs(title=paste(str_cell_type))
  ggsave(
    plot = volcano,
    filename = paste0(opt$in_dir,"/volcano_plots/",str_cell_type,".png")
    )

  # Top N DEG by p-value regardless of log2FC
  ls_topN_DEG_by_pval <- df_res %>%
    dplyr::filter(padj < .05) %>%
    dplyr::arrange(padj) %>%
    head(topN) %>%
    dplyr::pull(X)

  # Top N DEG ranked by positive log2FC
  ls_topN_DEG_up_reg <- df_res %>%
    dplyr::filter(padj < .05) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    head(topN) %>%
    dplyr::pull(X)

  # Top N DEG ranked by negative log2FC
  ls_topN_DEG_down_reg <- df_res %>%
    dplyr::filter(padj < .05) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    tail(topN) %>%
    dplyr::pull(X)

  return(list("ls_topN_DEG_by_pval"= ls_topN_DEG_by_pval,
              "ls_topN_DEG_up_reg" = ls_topN_DEG_up_reg,
              "ls_topN_DEG_down_reg" = ls_topN_DEG_down_reg))

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

list2heatmap <- function(ls_genes, title, output_name ){
  df_DEG_log2tpm_batch_corrrected <- df_log2tpm_batch_corrrected %>%
    dplyr::filter(X %in% ls_genes) %>%
    tibble::column_to_rownames(var  = "X")

  df_sample_annotations <- df_metadata %>%
  # dplyr::select(Run,sigil_general, LM22) %>%
    dplyr::select(Run, LM22, data_source) %>%
    tibble::column_to_rownames("Run")

  up_deg_heatmap <- pheatmap(
    main = paste0(title),
    df_DEG_log2tpm_batch_corrrected,
    scale = "row",
    show_rownames=F,
    show_colnames=F,
    annotation_col=df_sample_annotations)

  save_pheatmap_pdf(
    up_deg_heatmap,
    paste0(opt$out_dir,output_name,".pdf"))

}

make_umap <- function(num_neighbor,meta_col,df,out_path) {

  set.seed(123)
  df <- df %>%
    dplyr::select(-X)

  # Drop genes with low variance.
  getVar <- apply(df[, -1], 1, var)
  param <- median(getVar)
  log2trans_dat_filt <- df[getVar > param & !is.na(getVar), ]

  # Transpose and format
  log2trans_dat_filt_t <- as.data.frame(t(log2trans_dat_filt)) %>%
    dplyr::filter(!row.names(.) %in% c("gene"))
  rownames(log2trans_dat_filt_t) <- colnames(log2trans_dat_filt)

  # PCA.
  prcomp.out = as.data.frame(prcomp(log2trans_dat_filt_t, scale = F)$x)
  prcomp.out$Run = rownames(log2trans_dat_filt_t)
  # prcomp.out.merge = merge(prcomp.out, y = df)

  # Make color palette
  n <- length(unique(df_metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, df_metadata)

  # Plot UMAP
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP Kallisto: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste(out_path,meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)

}
###################
# MAIN
###################
# Arguments
option_list <- list(

  optparse::make_option(
    c("-i", "--in_dir"),
    type = "character",
    default = NULL,
    help = "path to inputs from combineGeneDE.R script"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "path to write outputs")
  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/refMatrix/"))){
  dir.create(file.path(opt$out_dir,"/refMatrix/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/volcano_plots/"))){
  dir.create(file.path(opt$out_dir,"/volcano_plots/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/UMAPs_DEG/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs_DEG/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/heatmaps_DEG/"))){
  dir.create(file.path(opt$out_dir,"/heatmaps_DEG/"),
              recursive = TRUE, showWarnings = TRUE)}

# Open files
df_metadata <- read.csv(file = paste0(opt$in_dir,"/metadata.csv"), header=TRUE)
df_log2tpm_batch_corrrected <- read.csv(
  file = paste0(opt$in_dir,"/combined_kallisto_log2tpm_batch_corrected.csv"),
  header=TRUE)

# Locate the deseq2 outputs
ls_deseq_paths = list.files(
  path = paste0(opt$in_dir,"/deseq2_outputs"),
  pattern = "*" , full.names = TRUE)
stopifnot(length(unique(df_metadata$LM22))==length(ls_deseq_paths))

# Import deseq2 outputs, make volcano, and get list of top up DEG
ls_topN_results <-lapply(ls_deseq_paths, read_in_deseq2, topN=20)

# Split DEG lists from read_in_deseq2 into different lists
ls_res_topN_DEG_by_pval <- unlist(lapply(c(1:length(ls_topN_results)),
    function(x) ls_topN_results[[x]]$ls_topN_DEG_by_pval))
ls_res_topN_upDEG <- unlist(lapply(c(1:length(ls_topN_results)),
    function(x) ls_topN_results[[x]]$ls_topN_DEG_up_reg))
ls_res_topN_downDEG <- unlist(lapply(c(1:length(ls_topN_results)),
    function(x) ls_topN_results[[x]]$ls_topN_DEG_down_reg))

##################################
# UMAPs and heatmaps using DEG
##################################

# Using all genes for reference
lapply(c(25,30), make_umap, meta_col="LM22",
    df = df_log2tpm_batch_corrrected, out_path = "UMAPs_DEG/all_genes")

# Filter df to top DEG and UMAP and heatmap
df_log2tpm_batch_corrrected_topDEG <- df_log2tpm_batch_corrrected %>%
    dplyr::filter(X %in% ls_res_topN_DEG_by_pval)
lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22",
    df = df_log2tpm_batch_corrrected_topDEG, out_path = "UMAPs_DEG/top_DEG_by_pval")
lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source",
    df = df_log2tpm_batch_corrrected_topDEG, out_path = "UMAPs_DEG/top_DEG_by_pval")
list2heatmap(ls_res_topN_DEG_by_pval, "All top DEG","/heatmaps_DEG/top_DEG_by_pval" )

# Filter df to top upregulated DEG and UMAP and heatmap
df_log2tpm_batch_corrrected_top_up_DEG <- df_log2tpm_batch_corrrected %>%
    dplyr::filter(X %in% ls_res_topN_upDEG)
lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22",
    df = df_log2tpm_batch_corrrected_top_up_DEG, out_path = "UMAPs_DEG/top_up_DEG")
lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source",
    df = df_log2tpm_batch_corrrected_top_up_DEG, out_path = "UMAPs_DEG/top_up_DEG")
list2heatmap(ls_res_topN_upDEG, "All top upregulated DEG", "/heatmaps_DEG/top_up_DEG")

# Filter df to top downregulated DEG and UMAP and heatmap
df_log2tpm_batch_corrrected_top_down_DEG <- df_log2tpm_batch_corrrected %>%
    dplyr::filter(X %in% ls_res_topN_downDEG)
lapply(c(5,10,15,20,25,30,35), make_umap, meta_col="LM22",
    df = df_log2tpm_batch_corrrected_top_down_DEG, out_path = "UMAPs_DEG/top_down_DEG")
lapply(c(5,10,15,20,25,30,35), make_umap, meta_col="data_source",
    df = df_log2tpm_batch_corrrected_top_down_DEG, out_path = "UMAPs_DEG/top_down_DEG")
list2heatmap(ls_res_topN_downDEG, "All top downregulated DEG", "/heatmaps_DEG/top_down_DEG")
