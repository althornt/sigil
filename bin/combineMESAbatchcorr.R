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

# Open files
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(head(df_manifest))




#log2trans_dat <- as.data.frame(log2(dat +1))


# Log log transform


# Batch correct
