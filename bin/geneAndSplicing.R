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
library(UpSetR)

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

print(opt$spliceDir)
print(opt$geneDir)

# Read in files 

# Reference matrix 
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
    select( -column_label2,-group, -column_label) %>%
    rename(gene = overlapping)

print("lm6 splice z gene name")
print(head(df_lm6_splice_z_gene))
print(dim(df_lm6_splice_z_gene))

df_lm6_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE) %>%
    rename(gene = X) 
print("lm6 gene")
print(head(df_lm6_gene_z))
print(dim(df_lm6_gene_z))

# merge z_score dfs 
df_merged_z <- df_lm6_gene_z %>%
        inner_join(df_lm6_splice_z_gene, by = "gene", suffix=c("_splice","_gene"))
# print(head(df_lm6_gene_z))
# print(head(df_lm6_splice_z_gene))
print("merged")
print(head(df_merged_z))
print(dim(df_merged_z))

# df_merged_z_long <- df_merged_z %>%
#     tidyr::gather(., key="idk", value = z, -c(gene, cell_type)) %>%
#     as.data.frame()
# print("df_lm6_gene_z_long")
# print(head(df_merged_z_long))


#####################################
# Scatter plots gene vs splicing 
#####################################
# Formatting long dfs
df_lm6_gene_z_long <- df_lm6_gene_z %>%
    tidyr::gather(., key="cell_type", value = gene_z, -c(gene)) %>%
    as.data.frame()
  
print("df_lm6_gene_z_long")
print(head(df_lm6_gene_z_long))

df_lm6_splice_z_long <- df_lm6_splice_z_gene %>%
    rename(cell_type_sig_splice = cell_type) %>%
    tidyr::gather(., key="cell_type", value = splice_z, -c(gene, event, cell_type_sig_splice)) %>%
    as.data.frame()

print("--------------")
print("df_lm6_splice_z_long")
print(head(df_lm6_splice_z_long))

# Merging the long dfs 
df_merged_z_long <-  df_lm6_splice_z_long %>%
        inner_join(df_lm6_gene_z_long, by = c("gene","cell_type"), suffix=c("_splice","_gene")) %>% 
         tidyr::drop_na()
print("merged z long")
print(head(df_merged_z_long))
print(dim(df_merged_z_long))


# Scatter plot with all matching events/genes
p <- ggplot(aes(x=gene_z, y=splice_z, color = cell_type), data=df_merged_z_long)+ 
  geom_point() +
  geom_text(
            label= df_merged_z_long$gene,
            nudge_x = 0.05, nudge_y = 0.05,
            check_overlap =F, col = "black", size = 1.5
          )
ggsave(plot = p, filename = paste0(opt$out_dir, "/lm6_scatter_color_by_splice_cell_type_all.png"))


# Make scatter plot for each matching gene
ls_genes <- unique(df_merged_z_long$gene)
for (g in ls_genes){
  df_g <- df_merged_z_long %>%
    filter(gene == g)
  p <- ggplot(aes(x=gene_z, y=splice_z, color = cell_type), data=df_g)+ 
    geom_point() +
    labs(title= g) + geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
    geom_text(
              label= df_g$cell_type_sig_splice,
              nudge_x = 0.05, nudge_y = 0.05,
              check_overlap =F, col = "black", size = 1.5
            ) +
    theme_minimal()
  ggsave(plot = p, filename = paste0(opt$out_dir, "/lm6_scatter_color_by_splice_cell_type_", g, ".png"))
}


