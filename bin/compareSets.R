#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
# library(uwot)
library(RColorBrewer)
# library(foreach)
library(doParallel)
# cl <- makeCluster(detectCores() - 1, outfile = "")
# registerDoParallel(cl)
library(uwot)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(purrr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)

# Arguments
option_list <- list(
  optparse::make_option(
    c("-s", "--splice_set"),
    type = "character",
    default = NULL,
    help = ""),

  optparse::make_option(
    c("-d", "--splice_set_df"),
    type = "character",
    default = NULL,
    help = ""),

  optparse::make_option(
    c("-g", "--gene_set_df"),
    type = "character",
    default = NULL,
    help = ""),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full "),
  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs")
    )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
if (!dir.exists(paste0(opt$out_dir))){
  dir.create(paste0(opt$out_dir),
  recursive = TRUE, showWarnings = TRUE)
}

# Open files
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  arrange(Run) %>%
  dplyr::select(Run,main_label,group_label, sigil_general) %>%
  # dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Read in and sort metadata 
metadata <- metadata %>% 
    arrange(Run) %>%
    tibble::column_to_rownames("Run") %>%
    select(main_label, group_label)
    # select(main_label, group_label, data_source)

# Splice gene set
df_spice_set <- read.csv(file = opt$splice_set_df)

print(head(df_spice_set))

df_gene_set <- read.csv(file = opt$gene_set_df)
print(head(df_gene_set))

# Get difference between splice and gene set 
cat("\n")
cat("Sets in gene set but not splice set \n ")
droplevels(unique(df_gene_set$set)[!(unique(df_gene_set$set) %in% unique(df_spice_set$set))])
cat("\n")

cat("Sets in splice set but not gene set \n ")
droplevels(unique(df_spice_set$set)[!(unique(df_spice_set$set) %in% unique(df_gene_set$set))])
cat("\n")

# Get union of sets in both to iterate through 
ls_all_sets <- union(df_spice_set$set, df_gene_set$set)
# print(length(ls_all_sets))

# Get intersection of spliced genes and genes
gene_intersection <- intersect(unique(df_spice_set$overlapping), unique(df_gene_set$X))
# print(gene_intersection)

cat("\n Number of common genes in gene and splice sets: \n")
print(length(gene_intersection))

cat("\n Number of unique genes in splice set: \n")
print(length(unique(df_spice_set$overlapping)))

cat("\n Number of unique junctions in splice set: \n")
print(length(unique(df_spice_set$event)))

cat("\n Number of unique genes in gene set: \n")
print(length(unique(df_gene_set$X)))
cat("\n")

# Compare on set level 
df_counts <- data.frame(matrix(ncol = 2, nrow = length(ls_all_sets)))
rownames(df_counts) <- ls_all_sets
colnames(df_counts) <- list("SpliceCount", "GeneCount")

for (s in ls_all_sets){
  print(s)

  # Counting 
  cnt_splice <- nrow(df_spice_set %>% 
    filter(set == s) )
  cnt_gene <- nrow(df_gene_set %>% 
    filter(set == s) )
  df_counts[s, "SpliceCount"] <- cnt_splice
  df_counts[s, "GeneCount"] <- cnt_gene

  # Comparing 



}


print(df_counts)

# Loop through each set 

    # count set gene 

    # count set splice

    # find intersection of genes  




# also llop through compare combined UP and DN - drop numbers from X or _UP from set 