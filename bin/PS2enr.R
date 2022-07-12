#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(pheatmap)
library(purrr)
library(tidyr)
library(stringr)

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--inputPS"),
    type = "character",
    default = NULL,
    help = " "),

  optparse::make_option(
    c("-o", "--outputPS"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"), 

  optparse::make_option(
    c("-c", "--removeChr"),
    type = "character",
    action = "store_true",
    default = NULL,
    help = "full path to put outputs")
  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Read in all MESA PS 
df_all_PS <- read.table(file = opt$inputPS, stringsAsFactors=FALSE,
                          sep="\t", header = TRUE, row.names=1) 
# print(head(df_all_PS))
print(dim(df_all_PS))

# Drop if all rows are NaN
df_all_PS <- df_all_PS[rowSums(is.na(df_all_PS)) != ncol(df_all_PS), ] 
# print(head(df_all_PS))
print(dim(df_all_PS))

df_all_PS <- as.data.frame(df_all_PS *100)
df_all_PS <- as.data.frame(log2(df_all_PS +1))
# print(head(df_all_PS))
print(dim(df_all_PS))

if (opt$removeChr){
  #reformat events/rownames to drop "chr"
  rownames(df_all_PS) <- stringr::str_remove(rownames(df_all_PS) , "chr")
}


# Write output
write.table(
  df_all_PS,
  file.path(opt$outputPS),
  quote=F,sep="\t", na="nan", 
  col.names = NA, row.names= TRUE)