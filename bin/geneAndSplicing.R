#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(uwot)
library(RColorBrewer)
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

print("hello")

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--spliceDir"),
    type = "character",
    default = NULL,
    help = " "),

  optparse::make_option(
    c("-c", "--geneDir"),
    type = "character",
    default = NULL,
    help = " "),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# # Open files
# metadata <- read.csv(file = opt$metadata)
# df_sample_annotations <- metadata %>%
#   dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
#   tibble::column_to_rownames("Run") %>%
#   t()

print(opt$spliceDir)
print(opt$geneDir)

df_event2gene <- read.csv(file= paste0(opt$spliceDir,"/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) 
print("splice ref mat")
print(head(df_event2gene))
print(dim(df_event2gene))

df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE) %>%
      rename(event = X)
print("lm6 splice z")
print(head(df_lm6_splice_z))
print(dim(df_lm6_splice_z))

# Add gene name to splice events
df_lm6_splice_z_gene <- merge(x=df_event2gene,y=df_lm6_splice_z,by="event",all=TRUE) %>%
    select(-event, -column_label2,-group, -cell_type, -column_label) %>%
    rename(gene = overlapping)

print("lm6 splice z gene name")
print(head(df_lm6_splice_z_gene))
print(dim(df_lm6_splice_z_gene))
print(head(df_lm6_splice_z_gene$gene))


df_lm6_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE)
print(head(df_lm6_gene_z))

df_lm6_gene_z_long <- df_lm6_gene_z %>%
    tidyr::gather(., key="sample", value = gene_z, -c(X)) %>%
    rename(gene = X) %>%
    as.data.frame()


print(head(df_lm6_gene_z_long))

df_lm6_splice_z_long <- df_lm6_splice_z_gene %>%
    tidyr::gather(., key="sample", value = splice_z, -c(gene)) %>%
    as.data.frame()

print(head(df_lm6_splice_z_long))

print("----")

# df_merged_z <-  merge(x=df_lm6_gene_z_long,y=df_lm6_splice_z_long,by="gene",all=TRUE) %>%
df_merged_z <-  df_lm6_splice_z_long %>%
        inner_join(df_lm6_gene_z_long, by = "gene", suffix=c("_splice","_gene")) %>% 
         tidyr::drop_na()
print(head(df_merged_z))

p <- ggplot(aes(x=gene_z, y=splice_z, color = sample_splice), data=df_merged_z)+ geom_point() 
ggsave(plot = p, filename = paste0(opt$out_dir, "/lm6_scatter_col_by_splice.png"))

p <- ggplot(aes(x=gene_z, y=splice_z, color = sample_gene), data=df_merged_z)+ geom_point() 
ggsave(plot = p, filename = paste0(opt$out_dir, "/lm6_scatter_col_by_gene.png"))


df_merged_z_matched <- df_merged_z %>%
    filter(sample_splice == sample_gene)

p <- ggplot(aes(x=gene_z, y=splice_z, color = sample_gene), data=df_merged_z_matched)+ 
geom_point(size = 1) +
geom_text(
            label= df_merged_z_matched$gene,
            nudge_x = 0.05, nudge_y = 0.05,
            check_overlap =F, col = "black", size = 1.5
          )
ggsave(plot = p, filename = paste0(opt$out_dir, "/lm6_scatter_col_by_gene_matched_cell_type.png"))
