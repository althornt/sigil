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
library(UpSetR)


sigil_out_path = "/mnt/results_sigil_combine/sigil_results_SongChoi_newlabels_20220614/"
output_path = "/mnt/figures/"


# metadata = 
df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"))

print(df_metadata)

#####################################
# Fig1B + C cell type summary plot
######################################

metadata_summary <- df_metadata %>%
  dplyr::arrange(main_label) %>%
  dplyr::count(main_label, data_source) 

print(metadata_summary)

ggplot(df_metadata, aes(x=group_label, fill=data_source)) +
  geom_bar(position="stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="None",    axis.title.x = element_blank()) +
  scale_fill_manual(values=c("#999999", "black")) 

ggsave(paste0(output_path,"metadata_bar_group.png"), width=10, height=5, unit="cm")

ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
  geom_bar(position="stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="None",    axis.title.x = element_blank()) +
    scale_fill_manual(values=c("#999999", "black"))

ggsave(paste0(output_path,"metadata_bar_main.png"), width=10, height=10, unit="cm")

ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
  geom_bar(position="stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="top",    axis.title.x = element_blank()) +
    scale_fill_manual("data source", values=c("#999999", "black"))

ggsave(paste0(output_path,"metadata_bar_main_legend.png"), width=10, height=10, unit="cm")

#####################################
# Fig1D UpSet plot 
######################################

# Read in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(sigil_out_path,"/combine_gene_out/ref_matrix/combined_main_group_withingroup_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) %>%
              rename(gene = X)

# Read in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                               sigil_out_path,
                                "/combine_mesa_out/ref_matrix_PS/combined_main_group_withingroup_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 

print(head(df_splice_ref))
print(dim(df_splice_ref))

# Read in MESA intron retention reference matrix
df_IR_ref <- read.csv(file= paste0(
                                sigil_out_path,
                                "/combine_mesa_out/ref_matrix_IR/combined_main_group_withingroup_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 
print(head(df_IR_ref))
print(dim(df_IR_ref))

# Genes affected by splice, IR, and gene
ls_gene_ref_genes <- unique(df_gene_ref$gene)
ls_splice_ref_genes <- unique(df_splice_ref$overlapping)
ls_IR_ref_genes <- unique(df_IR_ref$overlapping)

# Find genes in all 3
ls_ref_genes_splice_IR <- Reduce(intersect, 
                          list(ls_gene_ref_genes,
                              ls_splice_ref_genes,
                              ls_IR_ref_genes))
print("Intersection of all:")
print(ls_ref_genes_splice_IR)
for (i in ls_ref_genes_splice_IR){

    cat("\n",i,"\n")
}

print("Number of gene ref genes:")
print(length(ls_gene_ref_genes))
print("Number of splice ref genes:")
print(length(ls_splice_ref_genes))
print("Number of IR ref genes:")
print(length(ls_IR_ref_genes))

listInput <- list("Gene" = ls_gene_ref_genes,
                 "Splice" = ls_splice_ref_genes, 
                 "Intron Retention" = ls_IR_ref_genes)

plotObject <- upset(fromList(listInput), order.by = "freq", text.scale = c(1.6, 1.6, 1.6, 1.6, 2, 1.3),sets.bar.color=c("orange","skyblue","darkgreen"))
png(file= paste0(output_path, "/upsetplot_ref_genes_vs_splice_vs_IR.png"),
    width = 20,
    height    = 14,
    units     = "cm",
    res       = 1200)
print(plotObject)
dev.off()


# up_mat <- ComplexHeatmap::make_comb_mat(ComplexHeatmap::list_to_matrix(listInput),
#                                                     mode = "union",
#                                                     remove_empty_comb_set = TRUE)
# # up_mat_2 <- ComplexHeatmap::normalize_comb_mat(up_mat, full_comb_sets = TRUE)

# print(up_mat)

# png(file=paste0(output_path,"/upset.png"),
#     width = 10,
#     height    = 15,
#     units     = "cm",
#     res       = 1200)

# # png(file=paste0(output_path,"/upset.png"))
# up <- ComplexHeatmap::UpSet(up_mat, top_annotation=upset_top_annotation(up_mat, add_numbers = T))

# # # up <- ComplexHeatmap::UpSet(up_mat, comb_order = order(comb_size(up_mat)))

# # print(up)

# draw(up)
# dev.off()

