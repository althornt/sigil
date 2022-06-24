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
# library(reshape2)
library(ComplexHeatmap)
library(UpSetR)
library(cowplot)


sigil_out_path = "/mnt/results_sigil_combine/sigil_results_SongChoi_newlabels_20220614/"
fig_output_path = "/mnt/figures/"
if (!dir.exists(fig_output_path)){
dir.create(fig_output_path,
recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"), stringsAsFactors=FALSE)
print(df_metadata)

df_splice_ref_z <- read.csv(file = paste0(sigil_out_path, "gene_and_splicing_out/splice_ref_zscores.csv"))

# Read in all MESA PS 
df_all_PS <- read.table(file = paste0(sigil_out_path,"combine_mesa_out/batch_corr_mesa_allPS.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)
# Read in all gene exp 
df_exp <- read.table(file= paste0(sigil_out_path,
                    "combine_gene_out/combined_kallisto_log2tpm_batch_corrected.csv"),
                    sep=",", header = TRUE, row.names=1)
df_exp <- df_exp %>% mutate_if(is.character,as.numeric)



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

# ggsave(paste0(fig_output_path,"metadata_bar_group.png"), width=10, height=5, unit="cm")

# ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.position="None",    axis.title.x = element_blank()) +
#     scale_fill_manual(values=c("#999999", "black"))

# ggsave(paste0(fig_output_path,"metadata_bar_main.png"), width=10, height=10, unit="cm")

# ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.position="top",    axis.title.x = element_blank()) +
#     scale_fill_manual("data source", values=c("#999999", "black"))

# ggsave(paste0(fig_output_path,"metadata_bar_main_legend.png"), width=10, height=10, unit="cm")

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
# png(file= paste0(fig_output_path, "/upsetplot_ref_genes_vs_splice_vs_IR.png"),
#     width = 20,
#     height    = 14,
#     units     = "cm",
#     res       = 1200)
# print(plotObject)
# dev.off()

##############################
# fig 2
#########################
if (!dir.exists(paste0(fig_output_path,"PS/"))){
dir.create(paste0(fig_output_path,"PS/"),
recursive = TRUE, showWarnings = TRUE)}

print(head(df_metadata))


plot_PS <- function(junc,gene, comp_type, cell_type ){
    samples_main <- samples_other <- NA
    print("____________________________________________")
    print(cell_type)
    # replace, not sure why needed twice
    cell_type_ <- sub("_", " ", cell_type)
    cell_type_ <- sub("_", " ", cell_type_)

    print(cell_type_)
    print(comp_type)
    print(junc)

    if (comp_type == "main_label"){
        samples_main <- df_metadata %>%
                    filter(main_label == cell_type_) %>%
                    pull(Run)
        samples_other <- df_metadata %>%
                    filter(main_label != cell_type_) %>%
                    pull(Run)

    } else if (comp_type == "group_label"){
        samples_main <- df_metadata %>%
                    filter(group_label == cell_type_) %>%
                    pull(Run)
        samples_other <- df_metadata %>%
                    filter(group_label != cell_type_) %>%
                    pull(Run)
    } else {
        # Within group comparison 

        #get group label for the given label  
        group <- df_metadata %>%
                filter(main_label == cell_type_) %>%
                pull(group_label)
    
        group <- unique(group)
        print(group)
        print("_")
        
        samples_main <- df_metadata %>%
                     filter(main_label == cell_type_) %>%
                     pull(Run)

        samples_other <- df_metadata %>%
                    filter(main_label != cell_type_) %>%
                    filter(group_label == eval(group)) %>%
                    pull(Run)
    }

    
    print(length(samples_main))
    print(length(samples_other))

    ##############
    # PS Plot
    #############

    df_PS_junc <- df_all_PS %>% 
        filter(row.names(.) %in% c(junc)) %>%
        t() %>%
        as.data.frame() 

    df_PS_junc_main <- df_PS_junc %>%
        filter(row.names(.) %in% samples_main) %>%
        mutate(type = cell_type)
    df_PS_junc_other <- df_PS_junc %>%
        filter(row.names(.) %in% samples_other) %>%
        mutate(type = "other")

    df_PS_junc_p <- rbind(df_PS_junc_main,df_PS_junc_other)
    names(df_PS_junc_p) <- c("junc","type")

    p1 <- ggplot(df_PS_junc_p, aes(x=type, y=junc)) +
        geom_boxplot(outlier.shape = NA) +
        theme_classic() +
        geom_point(alpha = .5,position = position_jitter(seed = 1, width = 0.2)) +
        # theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
        theme(axis.title.x = element_blank()) + 
        labs(y= "PS", title = junc) 

    ##############
    # Expression Plot
    #############
    df_gene <- df_exp %>% 
        filter(row.names(.) %in% c(gene)) %>%
        t() %>%
        as.data.frame() 

    df_gene_main <- df_gene %>%
        filter(row.names(.) %in% samples_main) %>%
        mutate(type = cell_type)
    df_gene_other <- df_gene %>%
        filter(row.names(.) %in% samples_other) %>%
        mutate(type = "other")

    df_gene_p <- rbind(df_gene_main,df_gene_other)
    names(df_gene_p) <- c("junc","type")
    p2 <- ggplot(df_gene_p,aes(x=type, y=junc)) +
        geom_boxplot(outlier.shape = NA) +
        theme_classic() +
        geom_point(alpha = .5, position = position_jitter(seed = 1, width = 0.2)) +
        # theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        theme(axis.title.x = element_blank()) + 
        labs(y= "log2(TPM)", title =gene) +
        scale_y_continuous(
            expand = c(0,0),
            limits = c(0,  max(df_gene_p$junc) )) +
        coord_cartesian(clip = "off") 


    plot_grid(p1, p2)
    ggsave(paste0(fig_output_path,"PS/",junc,"_",gene,"_",cell_type,"_",comp_type,"_genePS.png"), width=15, height=8, unit="cm")



}

# for (row in 1:nrow(df_splice_ref_z)) {
for (row in 1:30) {
    junc <- as.character(df_splice_ref_z[row, "event"])
    gene <- as.character(df_splice_ref_z[row, "overlapping"])
    comp_type <- as.character(df_splice_ref_z[row, "group"])
    cell_type <- as.character(df_splice_ref_z[row, "cell_type"])

    plot_PS(junc,gene, comp_type, cell_type )

}



# df_splice_ref_z %>%
#     arrange(desc(ratio)) %>%
#     head(100)