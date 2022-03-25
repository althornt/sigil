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


########################
# Functions
#########################
format_merge <- function(df_gene,  df_splice){
  # Making long dfs
  df_gene_z_long <- df_gene %>%
      tidyr::gather(., key="cell_ID", value = gene_z, -c(gene, cell_types)) %>%
      as.data.frame()

  df_splice_z_long <- df_splice %>%
      tidyr::gather(., key="cell_ID", value = splice_z, -c(gene, event, cell_types)) %>%
      as.data.frame()

  # Merging the long dfs 
  df_merged_z_long <-  df_splice_z_long %>%
          inner_join(df_gene_z_long, by = c("gene","cell_ID"), suffix=c("_splice","_gene")) %>% 
          tidyr::drop_na()

  df_merged_z_long$splice_z = as.numeric(as.character(df_merged_z_long$splice_z))
  df_merged_z_long$gene_z = as.numeric(as.character(df_merged_z_long$gene_z))

  return(df_merged_z_long)

}
scatter_per_gene <- function(df_long, ouput_prefix){

  # Make scatter plot for each matching gene
  ls_genes <- unique(df_long$gene)
  for (g in ls_genes){
    df_g <- df_long %>%
      filter(gene == g)
    p <- ggplot(aes(x=gene_z, y=splice_z, color = cell_ID, shape = event), data=df_g)+ 
      geom_point() +
      labs(title= g) + geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
      geom_text(
                label= df_g$cell_types_splice,
                nudge_x = 0.05, nudge_y = 0.05,
                check_overlap =F, col = "black", size = 1
              ) +
      theme_minimal()

    ggsave(plot = p, filename = paste0(ouput_prefix, g, ".png"))
  }

}


scatter_all <- function(df_long, output_path){
  # Scatter plot with all matching events/genes
  p <- ggplot(aes(x=gene_z, y=splice_z, color = gene), data=df_long)+ 
    geom_point() +
    geom_text(
              label= df_long$gene,
              nudge_x = 0.05, nudge_y = 0.05,
              check_overlap =F, col = "black", size = 1.5
            )
  ggsave(plot = p, filename = output_path)
}

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

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/LM6_scatterplots/"))){
  dir.create(file.path(opt$out_dir,"/LM6_scatterplots/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/LM22_scatterplots/"))){
  dir.create(file.path(opt$out_dir,"/LM22_scatterplots/"),
              recursive = TRUE, showWarnings = TRUE)}

# Read in files 
# Events in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 
print("splice ref mat")
print(head(df_splice_ref))
print(dim(df_splice_ref))

print(unique(df_splice_ref$group))

# Map events to their signifcant cell type
df_event2cell <- df_splice_ref %>%
    group_by(event) %>%
    summarize(context = list(cell_type)) %>%
    mutate(cell_types = map_chr(context, toString)) %>%
    select(event, cell_types) %>%
    as.data.frame()
print(head(df_event2cell))
print(dim(df_event2cell))

df_event2cell <- merge(x = df_event2cell, y = df_splice_ref, by="event") %>%
  distinct(event, .keep_all = TRUE)

print(head(df_event2cell))
print(dim(df_event2cell))

print("=====================================================================")

# Genes in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(
                              opt$geneDir,
                              "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) %>%
              rename(gene = X)
print("gene ref mat")
print(head(df_gene_ref))

print(unique(df_gene_ref$group))
print(dim(df_gene_ref))


# Map genes to to their signifcant cell type
df_DEG2cell <- df_gene_ref %>%
    group_by(gene) %>%
    summarize(context = list(cell_type)) %>%
    mutate(cell_types = map_chr(context, toString)) %>%
    select(gene, cell_types) %>%
    as.data.frame()

print(head(df_DEG2cell))
print(dim(df_DEG2cell))

#######################
# LM6 group z_scores
########################

# Open LM6 splice z scores 
df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix/LM6_med_zscore.csv"),
                            header=TRUE) %>%
                    rename(event = X)

# Add sig cell type and gene name
df_lm6_splice_z <-  merge(x=df_lm6_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
    select( -column_label2,-group, -column_label, -cell_type) %>%
    rename(gene = overlapping)

print("----------------------------------------")
print(head(df_lm6_splice_z))
print(dim(df_lm6_splice_z))

# Open LM6 gene z scores 
df_lm6_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE) %>%
    rename(gene = X)
# Add  cell_type to z scores
df_lm6_gene_z <- merge(x=df_DEG2cell,y=df_lm6_gene_z,by="gene",all=TRUE) 
print("df_lm6_gene_z")
print(head(df_lm6_gene_z))
print(dim(df_lm6_gene_z))


# # Scatterplots 
df_lm6_merged_z_long <- format_merge(df_lm6_gene_z,df_lm6_splice_z )
print(head(df_lm6_merged_z_long))
scatter_all(df_lm6_merged_z_long, paste0(opt$out_dir, "/lm6_scatter_all.png"))
scatter_per_gene(df_lm6_merged_z_long, paste0(opt$out_dir, "/LM6_scatterplots/"))


#######################
# LM22
########################
# Open LM22 splice z scores 
df_lm22_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix/LM22_med_zscore.csv"),
                            header=TRUE) %>%
                    rename(event = X)

# Add sig cell type and gene name
df_lm22_splice_z <-  merge(x=df_lm22_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
    select( -column_label2,-group, -column_label, -cell_type) %>%
    rename(gene = overlapping)

print("----------------------------------------")
print(head(df_lm22_splice_z))
print(dim(df_lm22_splice_z))

# Open LM22 gene z scores 
df_lm22_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM22_med_zscore.csv"), header=TRUE) %>%
    rename(gene = X)
# Add  cell_type to z scores
df_lm22_gene_z <- merge(x=df_DEG2cell,y=df_lm22_gene_z,by="gene",all=TRUE) 
print("df_lm22_gene_z")
print(head(df_lm22_gene_z))
print(dim(df_lm22_gene_z))

# # Scatterplots 
df_lm22_merged_z_long <- format_merge(df_lm22_gene_z,df_lm22_splice_z )
scatter_all(df_lm22_merged_z_long, paste0(opt$out_dir, "/lm22_scatter_all.png"))
scatter_per_gene(df_lm22_merged_z_long, paste0(opt$out_dir, "/LM22_scatterplots/"))



#######################
# Count splicing types 
########################