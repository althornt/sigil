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
library(purrr)

importMetaMESA <- function(row){
  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general, LM22)

  # Add metadata to column
  df_metadata$data_source <- row[1] # add name of data source
  df_metadata$type <- row[4] # add rna-seq type (paired vs single)

  # Get paths to MESA inclusion count files
  res_path <- file.path(row[2], "mesa_out", "mesa_inclusionCounts.tsv")

  return(list(
      "ls_kallsto_paths"=res_path,
      "metadata"=df_metadata,
      "ls_samples_run"=df_metadata$Run))
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
if (!dir.exists(file.path(opt$out_dir,"/UMAPs_pre_batch_correction/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs_pre_batch_correction/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/UMAPs_post_batch_correction/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs_post_batch_correction/"),
              recursive = TRUE, showWarnings = TRUE)}

# Open manifest
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(head(df_manifest))

# Import and combine source metadata files
ls_mesa_meta = apply(df_manifest, 1, importMetaMESA)

# Split into kallisto and metadata files for each data set
ls_mesa_files <- ls_meta <- ls_sample_names <- c()
for (item in ls_mesa_meta) {
     ls_mesa_files <- append(ls_mesa_files, item[1])
     ls_meta <- append(ls_meta, item[2])
     ls_sample_names <- append(ls_sample_names, item[3])
   }

# Merge metadata from each data source by rows
df_merged_metadata <- do.call("rbind", ls_meta)
rownames(df_merged_metadata) <- c()

# Remove samples without LM22 labels from metadata
df_merged_metadata_lm22 <- df_merged_metadata %>%
   dplyr::filter(LM22 != "")
write.csv(df_merged_metadata_lm22,
            file.path(file.path(opt$out_dir,"metadata.csv")),
            row.names = FALSE)

# List of samples with LM22 labels
ls_smpls_lm22 <- as.character(df_merged_metadata_lm22$Run)
print(typeof(unlist(ls_smpls_lm22)))
print(unlist(ls_smpls_lm22))

# Merge mesa inclusion count files
df_mesa_merge <- unlist(ls_mesa_files) %>%
  lapply(read.csv, sep = "\t") %>%
  purrr::reduce(inner_join, by = "cluster")

# Drop non LM22 samples from mesa counts
df_mesa_merge_lm22 <- df_mesa_merge %>%
  dplyr::select(ls_smpls_lm22)

print(head(df_mesa_merge_lm22))
print(dim(df_mesa_merge_lm22))

#  Log2 + 1 transform
df_mesa_merge_lm22_log2 <- as.data.frame(log2(df_mesa_merge_lm22 +1))

head(df_mesa_merge_lm22_log2)
print(dim(df_mesa_merge_lm22_log2))

# Batch correction
df_mesa_merge_lm22_log2_batch_corr <- limma::removeBatchEffect(
                                  df_mesa_merge_lm22_log2,
                                  batch = df_merged_metadata_lm22$data_source,
                                  batch2 = df_merged_metadata_lm22$type
                                  )

head(df_mesa_merge_lm22_log2_batch_corr)
dim(df_mesa_merge_lm22_log2_batch_corr)
