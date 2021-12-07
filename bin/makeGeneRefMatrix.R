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

  # Get cell type name from file path
  str_cell_type <- tools::file_path_sans_ext(basename(file_path))

  df_res <- read.csv(file = file_path, header = TRUE)

  # Get DEG
  df_res_gene <- df_res %>%
    dplyr::filter(padj < .05) %>%
    dplyr::arrange(desc(log2FoldChange))

  # Get top 50 upregulated DEG
  ls_topN_upreg_DEG <- df_res_gene %>%
    head(topN) %>%
    dplyr::pull(X)

  # Make volcano plot
  volcano <- ggplot(data=df_res, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point() + theme_minimal() +
    geom_hline(yintercept=-log10(0.05), col="red") +
    labs(title=paste(str_cell_type))

  ggsave(plot = volcano,filename = paste0(opt$in_dir,"/volcano_plots/",str_cell_type,".png"))

  # return(list(eval(str_cell_type) = ls_top50_upreg_DEG))

  return(ls_topN_upreg_DEG)

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
ls_cell_topDEG <- lapply(ls_deseq_paths, read_in_deseq2, topN=10)
print(ls_cell_topDEG)

print("DEG overlap:")
print(Reduce(intersect, ls_cell_topDEG))

# Make heatmap using top up DEG from all cell types
ls_all_topDEG <-  unlist(ls_cell_topDEG, recursive = FALSE)

df_DEG_log2tpm_batch_corrrected <- df_log2tpm_batch_corrrected %>%
  dplyr::filter(X %in% ls_all_topDEG) %>%
  tibble::column_to_rownames(var  = "X")

df_sample_annotations <- df_metadata %>%
# dplyr::select(Run,sigil_general, LM22) %>%
  dplyr::select(Run, LM22) %>%
  tibble::column_to_rownames("Run")

up_deg_heatmap <- pheatmap(
  main = paste0(" All UP DEG "),
  df_DEG_log2tpm_batch_corrrected,
  scale = "row",
  show_rownames=F,
  show_colnames=F,
  annotation_col=df_sample_annotations)

save_pheatmap_pdf(
  up_deg_heatmap,
  paste0(opt$out_dir,"/UP_DEG_heatmap.pdf"))
