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
if (!dir.exists(output_path)){
dir.create(output_path,
recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"), stringsAsFactors=FALSE)
print(df_metadata)

df_splice_ref_z <- read.csv(file = paste0(sigil_out_path, "gene_and_splicing_out/splice_ref_zscores.csv"))

# Read in all MESA PS 
df_all_PS <- read.table(file = paste0(sigil_out_path,"combine_mesa_out/batch_corr_mesa_allPS.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)


# #####################################
# # Fig1B + C cell type summary plot
# ######################################

# metadata_summary <- df_metadata %>%
#   dplyr::arrange(main_label) %>%
#   dplyr::count(main_label, data_source) 

# print(metadata_summary)

# ggplot(df_metadata, aes(x=group_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(legend.position="None",    axis.title.x = element_blank()) +
#   scale_fill_manual(values=c("#999999", "black")) 

# ggsave(paste0(output_path,"metadata_bar_group.png"), width=10, height=5, unit="cm")

# ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.position="None",    axis.title.x = element_blank()) +
#     scale_fill_manual(values=c("#999999", "black"))

# ggsave(paste0(output_path,"metadata_bar_main.png"), width=10, height=10, unit="cm")

# ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.position="top",    axis.title.x = element_blank()) +
#     scale_fill_manual("data source", values=c("#999999", "black"))

# ggsave(paste0(output_path,"metadata_bar_main_legend.png"), width=10, height=10, unit="cm")

# #####################################
# # Fig1D UpSet plot 
# ######################################

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
# for (i in ls_ref_genes_splice_IR){
#     cat("\n",i,"\n")
# }

df_splice_ref %>%
    filter(overlapping %in% ls_ref_genes_splice_IR) %>%
    select(event, overlapping, cell_type, group)


df_IR_ref %>%
    filter(overlapping %in% ls_ref_genes_splice_IR) %>%
    select(event, overlapping, cell_type, group)

# for (i in ls_IR_ref_genes){
#     cat("\n",i)
# }


print("Number of gene ref genes:")
print(length(ls_gene_ref_genes))
print("Number of splice ref genes:")
print(length(ls_splice_ref_genes))
print("Number of IR ref genes:")
print(length(ls_IR_ref_genes))

# listInput <- list("Gene" = ls_gene_ref_genes,
#                  "Splice" = ls_splice_ref_genes, 
#                  "Intron Retention" = ls_IR_ref_genes)

# plotObject <- upset(fromList(listInput), order.by = "freq", text.scale = c(1.6, 1.6, 1.6, 1.6, 2, 1.3),sets.bar.color=c("orange","skyblue","darkgreen"))
# png(file= paste0(output_path, "/upsetplot_ref_genes_vs_splice_vs_IR.png"),
#     width = 20,
#     height    = 14,
#     units     = "cm",
#     res       = 1200)
# print(plotObject)
# dev.off()

##############################
# fig 2
#########################
if (!dir.exists(paste0(output_path,"PS/"))){
dir.create(paste0(output_path,"PS/"),
recursive = TRUE, showWarnings = TRUE)}



plot_PS <- function(junc, comp_type, cell_type ){

    # Get samples of main type
    samples_main <- df_metadata[df_metadata$main_label == "Monocyte NT",]$Run
    samples_other <- df_metadata[df_metadata$main_label != "Monocyte NT",]$Run 

    print(samples_main)
    print(samples_other)

    df_PS_junc <- df_all_PS %>% 
        filter(row.names(.) %in% c(junc)) %>%
        t() %>%
        as.data.frame() 

    df_PS_junc <- df_PS_junc %>%
        mutate(cell_type_match = ifelse(row.names(.) %in% samples_main, "main", "other"))

    ggplot(df_PS_junc, aes(x=cell_type_match, y=V1)) +
        geom_boxplot(outlier.shape = NA) +
        theme_classic() +
        geom_point(position = position_jitter(seed = 1, width = 0.2)) +
        theme(axis.title.x = element_blank()) + labs(y= "PS") 

    ggsave(paste0(output_path,"PS/",junc,".png"), width=10, height=10, unit="cm")


}

plot_PS("20:63570868-63572107:-", main_label, "Monocyte_NT"  )
plot_PS("12:57107314-57108162:-", main_label, "Monocyte_R848-18h"  )


# print(head(df_metadata))


df_splice_ref_z %>%
    arrange(desc(ratio)) %>%
    head(10)

# print(head(df_all_PS))
# print(head(df_metadata))