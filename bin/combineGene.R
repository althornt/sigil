#!/usr/bin/env Rscript

library(optparse)
library(magrittr)
library(dplyr)

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
    help = "path to tsv manifest file with 3 columns: data_source, res_path, metadata_path")

  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directory
if (!dir.exists(paste0(opt$out_dir,"/UMAPs/"))){
  dir.create(file.path(opt$out_dir), recursive = TRUE, showWarnings = FALSE)}

# Open files
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(head(df_manifest))

# Check all combined kallisto results files exist
ls_kallisto_paths <- file.path(
  df_manifest$res_path,
  "post_kallisto_out",
   "combined_kallisto_log2tpm.csv")
if(all(file.exists(ls_kallisto_paths)) != TRUE)
  stop("Error: missing combined_kallisto_log2tpm.csv files.")

importKallistoMeta <- function(row){
  source_name <- row[1]

  # Read kallisto file
  res_path <- file.path(row[2], "post_kallisto_out","combined_kallisto_log2tpm.csv")
  df_res <- read.csv(res_path) %>%
    dplyr::rename("gene" = "X")

  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general)

  df_metadata$data_source <- source_name
  print(df_metadata)

  return(list("res"=df_res, "metadata"=df_metadata))
}

# Import and combine kallisto files
ls_kallisto_meta = apply(df_manifest, 1, importKallistoMeta)

# Split into kallisto and metadata files for each data set
ls_kallisto <- ls_mesa <- c()
for (item in ls_kallisto_meta) {
     ls_kallisto <- append(ls_kallisto, item[1])
     ls_mesa <- append(ls_mesa, item[2])
   }

# Merge kallisto by columns
df_merged_kallisto_res <- Reduce(function(...) merge(..., all = TRUE, by="gene"), ls_kallisto)

# Merge metadata by rows
df_merged_metadata <- do.call("rbind", ls_mesa)
rownames(df_merged_metadata) <- c()
