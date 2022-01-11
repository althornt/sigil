#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--mesa_PS"),
    type = "character",
    default = NULL,
    help = " `mesa_allPS.tsv` input file "),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Open files
metadata = read.csv(file = opt$metadata)
all_PS = read.table(file = opt$mesa_PS, sep="\t", row.names = 1, header = TRUE)

# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/explore/"))){
  dir.create(paste0(opt$out_dir,"/explore/"),
   recursive = TRUE, showWarnings = TRUE)
}

print("test")

# get list of cell types
# read in file for each cell type
# heatmap top n events up and n down
# filter nan

# start with dendritic activated , resting
#be super strinngnne
#
