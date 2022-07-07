#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(uwot)
library(magrittr)
library(purrr)
library(tidyr)
library(ComplexHeatmap)
library(UpSetR)
library(cowplot)
# library(enrichR)
# library(GenomicRanges)
library(valr)
library(rGREAT)

sigil_out_path = "/mnt/results_sigil_combine/sigil_results_SongChoi_newlabels_20220614/"
fig_output_path = "/mnt/figures/"
if (!dir.exists(fig_output_path)){
    dir.create(fig_output_path,
    recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"), stringsAsFactors=FALSE)
# print(df_metadata)

df_splice_ref_z <- read.csv(file = paste0(sigil_out_path, "gene_and_splicing_out/splice_ref_zscores.csv"))
df_IR_ref_z <- read.csv(file = paste0(sigil_out_path, "gene_and_splicing_out/IR_ref_zscores.csv"))

# Read in all MESA PS 
df_all_PS <- read.table(file = paste0(sigil_out_path,"combine_mesa_out/batch_corr_mesa_allPS.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)
# Read in all gene exp 
df_exp <- read.table(file= paste0(sigil_out_path,
                    "combine_gene_out/combined_kallisto_log2tpm_batch_corrected.csv"),
                    sep=",", header = TRUE, row.names=1)
df_exp <- df_exp %>% mutate_if(is.character,as.numeric)

# Read in MESA intron retention 
df_IR_table <- read.table(file =paste0(sigil_out_path,
                          "/combine_mesa_out/batch_corr_mesa_ir_table_intron_retention.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_IR_table <- df_IR_table %>% mutate_if(is.character,as.numeric)

# Read in sigil sets
df_splice_set <- read.csv(file = paste0(sigil_out_path,
                          "combine_mesa_out/splice_set_PS/df_splice_sets.csv"))
print(head(df_splice_set))
df_IR_set <- read.csv(file = paste0(sigil_out_path,
                          "combine_mesa_out/splice_set_IR/df_splice_sets.csv"))
print(head(df_IR_set))

df_gene_set <- read.csv(file = paste0(sigil_out_path,
                          "combine_gene_out/gene_sets/df_gene_sets.csv"))
print(head(df_gene_set))


# df_bed_PS <- distinct(read.table(file = paste0(sigil_out_path,
#                         "combine_mesa_out/explore_ref_matrix_PS/RefMatrix.bed"), 
#                         header = FALSE)) 
# df_bed_IR <- distinct(read.table(file = paste0(sigil_out_path,
#                         "combine_mesa_out/explore_ref_matrix_IR/RefMatrix.bed"),
#                          header = FALSE))
# names(df_bed_PS) <-  names(df_bed_IR) <- list("CHROM", "START","STOP")                        

# UCSC alt event track table
df_alt_events <- read.table(file = "/mnt/files/UCSC_alt_events_track.tsv", header = TRUE)  %>%
    dplyr::rename(start = chromStart, end = chromEnd ) %>%
    mutate_at('start', as.integer) %>%
    mutate_at('end', as.integer)

print(head(df_alt_events))

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

# ###########################################
# # Fig1D UpSet plot comparing type of sets
# ############################################

# Read in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(sigil_out_path,"/combine_gene_out/ref_matrix/combined_main_group_withingroup_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) 
                         
                        #  %>%
            #   rename(gene = X)

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
ls_gene_ref_genes <- unique(df_gene_ref$X)
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

##########################################################
# Fig 2 compare junctions to UCSC Alt track alt promoter
############################################################
# Add bed cols to splice ref 
df_splice_set <- df_splice_set %>%
  rowwise() %>%
    mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1])) %>%
    mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
    mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
    mutate(strand = strsplit(event, ":")[[1]][3]) %>%
    as.data.frame() 
    # %>%
    # head(10)
print(head(df_splice_set))
print(dim(df_splice_set))
print("~~~~~~~~~~~")

# Get number to randomly sample from PS file 
num_unique_ref_juncs <- length(unique(df_splice_set$event))
print(length(unique(num_unique_ref_juncs)))

# Randonmly sample from the file with all junctions 
ls_random_sample <- sample(rownames(df_all_PS), num_unique_ref_juncs, replace=F)

print(head(ls_random_sample))
print(length(ls_random_sample))
print(length(unique(ls_random_sample)))

print("#############")

# Format random junctions into bed 
df_random_juncs <- as.data.frame(ls_random_sample) %>%
    dplyr::rename(event = ls_random_sample )  %>%
    rowwise() %>%
    mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1]))  %>%
    mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
    mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
    mutate(strand = strsplit(event, ":")[[1]][3]) %>%
    as.data.frame() 


print(head(df_random_juncs))
print(dim(df_random_juncs))


# Format all junctions into bed 
ls_all_juncs <- rownames(df_all_PS)
print(length(ls_all_juncs))
ls_all_juncs <- ls_all_juncs[!(ls_all_juncs %in% df_splice_set$event)]
print(length(ls_all_juncs))

df_all_background_juncs <- as.data.frame(ls_all_juncs) %>%
    dplyr::rename(event = ls_all_juncs )  %>%
    rowwise() %>%
    mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1]))  %>%
    mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
    mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
    mutate(strand = strsplit(event, ":")[[1]][3]) %>%
    as.data.frame() 


print(head(df_all_background_juncs))
print(dim(df_all_background_juncs))


# Set up alt promoter df
df_alt_events_promoter <- df_alt_events %>%
    filter(name=="altPromoter") %>%
    select("chrom","start","end","name","strand")


compAltProm <- function(df_bed ){

    # Find intersection of splice ref junctions with alt promoters; 
    # if they intersect/overlap it suggests the alt promoter its intersected with 
    # is being skipped
    df_alt_promoter_splice_ref_int <- bed_intersect(df_bed, df_alt_events_promoter, 
            suffix = c(".splice_ref", ".alt_pro"))  %>%
            as.data.frame() 
    print(head(df_alt_promoter_splice_ref_int))

    # Find which junctions within 200 of an alt promoter 
    # but do not overlap with an alt promoter, since thats accounted for in the other value
    df_alt_promoter_splice_ref_closest <- bed_closest(df_bed, df_alt_events_promoter, 
            suffix = c(".splice_ref", ".alt_pro"))  %>%
            as.data.frame() %>%
            filter(abs(.dist) <= 200) %>%
            filter(.overlap == 0) %>%
            as.data.frame()  
    head(df_alt_promoter_splice_ref_closest, n=20)

    # Count results
    unique_splice_juncs <- unique(df_bed$event)
    unique_splice_juncs_alt_pro_int <- unique(df_alt_promoter_splice_ref_int$event.splice_ref)
    unique_splice_juncs_alt_pro_close <- unique(df_alt_promoter_splice_ref_closest$event.splice_ref)
    unique_splice_juncs_alt_pro_either <- union(df_alt_promoter_splice_ref_int$event.splice_ref, 
                                                df_alt_promoter_splice_ref_closest$event.splice_ref )

    cat("\n Percent of splice ref junctions that contain an alt promoter: \n")
    perc_contain_alt <- length(unique_splice_juncs_alt_pro_int)/length(unique_splice_juncs)*100 
    cat(perc_contain_alt)
    cat("\n")

    cat("\n Percent of splice ref junctions that are 200 bp downstream of an alt promoter: \n")
    perc_downstream_alt <- length(unique_splice_juncs_alt_pro_close)/length(unique_splice_juncs)*100
    cat(perc_downstream_alt )
    cat("\n")

    cat("\n Percent of splice ref junctions that are either: \n")
    percent_alt <-length(unique_splice_juncs_alt_pro_either)/length(unique_splice_juncs)*100 
    cat(percent_alt)
    cat("\n")

    # Make alt pro calc res into df
    df_res_perc <- t(data.frame(list("Contains Alt Promoter"=perc_contain_alt,
                                    "Downstream Alt Promoter"=perc_downstream_alt, 
                                    "Alt Promoter"=percent_alt )))
    print(df_res_perc)
    df_res_counts <- t(data.frame(list("Contains Alt Promoter"= length(unique_splice_juncs_alt_pro_int),
                                    "Downstream Alt Promoter"=length(unique_splice_juncs_alt_pro_close),
                                     "Alt Promoter"=length(unique_splice_juncs_alt_pro_either) )))
    print(df_res_counts)
    # Add alt promoter info the input df 
    df_bed <- df_bed %>%
        mutate(contains_alt_pro = ifelse((event %in% unique_splice_juncs_alt_pro_int), 1, 0)) %>%
        mutate(downstream_alt_pro = ifelse((event %in% unique_splice_juncs_alt_pro_close), 1, 0)) %>%
        mutate(both = ifelse(((event %in% unique_splice_juncs_alt_pro_int) & (event %in% unique_splice_juncs_alt_pro_close)), 1, 0)) %>%
        mutate(either = ifelse(((event %in% unique_splice_juncs_alt_pro_int) | (event %in% unique_splice_juncs_alt_pro_close)), 1, 0)) 

    return(list(df_bed, df_res_perc, df_res_counts))

}

ls_alt_res_splice_sets <- compAltProm(df_splice_set)
df_splice_set <- ls_alt_res_splice_sets[[1]] 
head(df_splice_set)

print("~~~~~~~~~~~")
ls_alt_res_random <- compAltProm(df_random_juncs)
df_random_juncs <- ls_alt_res_random[[1]]
head(df_random_juncs)


print("~~~~~~~~~~~")
ls_alt_res_background <- compAltProm(df_all_background_juncs)
df_all_background <- ls_alt_res_background[[1]]
head(df_all_background)

# Make plot comparing alt promoter count in splice sets to random junctions

# df_alt_pro_res <- cbind(ls_alt_res_splice_sets[[2]], ls_alt_res_random[[2]], ls_alt_res_background[[2]])
# colnames(df_alt_pro_res) <- list("splice_ref", "random_sample", "background")


# Fishers exact test 

##                      Non-AFE      AFE
## Splice ref           4.5           4.5
## Non splice ref       2.5             2.5
df_alt_pro_res_count <- cbind(ls_alt_res_splice_sets[[3]], ls_alt_res_background[[3]])
colnames(df_alt_pro_res_count) <- list("splice_ref", "background")
print(df_alt_pro_res_count)


# Building contingency table 
AltPro_splice_ref <- df_alt_pro_res_count["Alt.Promoter","splice_ref"]
AltPro_splice_background <- df_alt_pro_res_count["Alt.Promoter","background"]

NotAltPro_splice_ref <- num_unique_ref_juncs - AltPro_splice_ref
NotAltPro_splice_background <- length(ls_all_juncs) - AltPro_splice_background

print(paste0(AltPro_splice_ref, AltPro_splice_background, NotAltPro_splice_ref, NotAltPro_splice_background, sep= " "))

df_fishers <- data.frame("nonAFE" = c(NotAltPro_splice_ref, NotAltPro_splice_background), 
                "AFE" = c(AltPro_splice_ref, AltPro_splice_background), 
                row.names = c("spliceRef", "background"))


# Run Fishers test
print(df_fishers)

res_fishers <- stats::fisher.test(df_fishers)
print(res_fishers)
print(res_fishers$p.value)

df_alt_pro_res_perc <- cbind(ls_alt_res_splice_sets[[2]], ls_alt_res_background[[2]])
colnames(df_alt_pro_res_perc) <- list("splice_ref", "background")
print(df_alt_pro_res_perc)


df_alt_pro_res_perc <- df_alt_pro_res_perc %>% 
    as.data.frame() %>%
    rownames_to_column( var= "type")
print(df_alt_pro_res_perc)

df_alt_pro_res_perc <- pivot_longer(df_alt_pro_res_perc,  values_to = "percent",  names_to = "name",  cols = c("splice_ref", "background")) 
print(df_alt_pro_res_perc)

# Number of cars in each class:
df_alt_pro_res_perc %>% 
    ggplot(aes(x=type, y=percent, fill=name)) +
    geom_bar(stat = 'identity', position = "dodge") +
    theme_classic() +
    scale_fill_manual("name", values=c("#999999", "black"))+
    theme(legend.title= element_blank(), axis.title.x=element_blank()) +
    labs(title = paste0("Fishers test p-value = ", signif(res_fishers$p.value, 4)))

ggsave(paste0(fig_output_path,"alt_promoter.png"), width=20, height=10, unit="cm", dpi = 300)



quit()


##########################################################
# Fig 2 rGREAT splice ref
############################################################

df_splice_set_clean <- df_splice_set %>%
    filter(! chrom %in% c("chrGL000219.1", "chrKI270711.1", "chrKI270745.1"))

# Run Great
job = submitGreatJob(makeGRangesFromDataFrame(df_splice_set_clean),
                        species= "hg38", version = "4")

# Get tables from Run
tb = getEnrichmentTables(job, category = c("GO"))


# print(head(tb[["GO Molecular Function"]], n=25))
# print(head(tb[[ "GO Biological Process" ]], n=50))

ls_great_plots <- list()
# for (pval_type in c("Hyper_Adjp_BH", "Binom_Adjp_BH","Binom_Fold_Enrichment" )) {
for (pval_type in c( "Binom_Adjp_BH" )) {

    # Make barplot of adjusted pvalues
    BP <- tb[["GO Biological Process"]] %>%
        mutate(pval_type = -log10(get(pval_type))) %>%
        arrange(desc(pval_type)) %>%
        head(40) %>%
        ggplot(aes(x = pval_type, y = reorder(name, pval_type)))  + 
            geom_bar(stat = 'identity') +
            theme_classic() +
            theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=20)) +
            labs(title = paste0("PS GO Biological Process"), x = "Binom_Adjp_BH")

    ggsave(paste0(fig_output_path,"rGREAT_GO_BP_all_splice_ref_",pval_type,".png"), width=30, height=20, unit="cm", dpi = 300)


    # Make barplot of adjusted pvalues
    MF <- tb[["GO Molecular Function"]] %>%
        mutate(pval_type = -log10(get(pval_type))) %>%
        arrange(desc(pval_type)) %>%
        head(40) %>%
        ggplot(aes(x = pval_type, y = reorder(name, pval_type)))  + 
             geom_bar(stat = 'identity') +
             theme_classic() +
             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=20)) +
             labs(title = paste0("PS GO Molecular Function"), x = "Binom_Adjp_BH")

    ggsave(paste0(fig_output_path,"rGREAT_GO_MF_all_splice_ref_", pval_type, ".png"), width=30, height=20, unit="cm", dpi = 300)

    plot_grid(BP, MF, nrow=1, align = 'v', axis = 'l')
    ggsave(paste0(fig_output_path,"rGREAT_GO_MF_BP_all_splice_ref_", pval_type, ".png"), width=70, height=30, unit="cm")

    ls_great_plots$IR_BP <- BP
    ls_great_plots$IR_MF <- MF

}

# # This isnt saving?
# plotObject <- plotRegionGeneAssociationGraphs(job)
# png(file= paste0(fig_output_path, "rGREAT_region_associations.png"),
#     width = 20,
#     height    = 14,
#     units     = "cm",
#     res       = 1200)
# print(plotObject)
# dev.off()

##########################################################
# Fig 2 rGREAT IR
############################################################
# df_IR_set <- df_IR_set %>%
#   rowwise() %>%
#     mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1])) %>%
#     mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
#     mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
#     mutate(strand = strsplit(event, ":")[[1]][3]) %>%
#     as.data.frame() 

# print(head(df_IR_set))


# # Run Great
# job_IR = submitGreatJob(makeGRangesFromDataFrame(df_IR_set),
#                         species= "hg38", version = "4")

# # Get tables from Run
# tb_IR = getEnrichmentTables(job_IR, category = c("GO"))


# print(head(tb_IR[["GO Molecular Function"]], n=25))
# print(head(tb_IR[[ "GO Biological Process" ]], n=60))

# for (pval_type in c("Binom_Adjp_BH" )) {
# # for (pval_type in c("Hyper_Adjp_BH", "Binom_Adjp_BH" )) {

#     # Make barplot of adjusted pvalues
#     BP <- tb_IR[["GO Biological Process"]] %>%
#         mutate(pval_type = -log10(get(pval_type))) %>%
#         arrange(desc(pval_type)) %>%
#         head(40) %>%
#         ggplot(aes(x = pval_type, y = reorder(name, pval_type)))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=20)) +
#             labs(title = paste0("IR GO Biological Process"), x = "Binom_Adjp_BH")
#     ggsave(paste0(fig_output_path,"rGREAT_GO_BP_all_IR_ref_",pval_type,".png"), width=30, height=20, unit="cm", dpi = 300)


#     # Make barplot of adjusted pvalues
#     MF <- tb_IR[["GO Molecular Function"]] %>%
#         mutate(pval_type = -log10(get(pval_type))) %>%
#         arrange(desc(pval_type)) %>%
#         head(40) %>%
#         ggplot(aes(x = pval_type, y = reorder(name, pval_type)))  + 
#              geom_bar(stat = 'identity') +
#              theme_classic() +
#              theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=20)) +
#              labs(title = paste0("IR GO Molecular Function"), x = "Binom_Adjp_BH")

#     ggsave(paste0(fig_output_path,"rGREAT_GO_MF_all_IR_ref_", pval_type, ".png"), width=30, height=20, unit="cm", dpi = 300)

#     plot_grid(BP, MF, nrow=1, align = 'v', axis = 'l')
#     ggsave(paste0(fig_output_path,"rGREAT_GO_MF_BP_all_IR_ref_", pval_type, ".png"), width=70, height=30, unit="cm")

#     ls_great_plots$PS_BP <- BP
#     ls_great_plots$PS_MF <- MF

# }

# print(length(ls_great_plots))
# plot_grid(ls_great_plots[[1]], ls_great_plots[[2]], ls_great_plots[[3]], ls_great_plots[[4]],
#              nrow=2, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"rGREAT_GO_MF_BP_all_PS_IR_ref_Binom_Adjp_BH.png"), width=75, height=60, unit="cm")


##########################################################
# Fig 2 Enrichr on all genes from each set type
############################################################
# if (!dir.exists(paste0(fig_output_path,"enrichr"))){
# dir.create(paste0(fig_output_path,"enrichr"),
# recursive = TRUE, showWarnings = TRUE)}


# dbs <- c("GO_Biological_Process_2021",
#         "GO_Cellular_Component_2021", 
#         "GO_Molecular_Function_2021",
#         "KEGG_2021_Human")

# ls_sigil_types <- list("gene"=ls_gene_ref_genes,
#                     "splice"=ls_splice_ref_genes,
#                     "IR"=ls_IR_ref_genes)
# ls_BP <-ls_CC <- ls_MF <- ls_kegg <- list()
# for (i in names(ls_sigil_types)){
#     cat("-------------------------------------------------------")
#     cat("\n",i,"\n")

#     # Run enrichr
#     enriched <- enrichr(ls_sigil_types[[i]], dbs)

#     # Save results to files
#     write.csv(enriched[[1]], paste0(fig_output_path,"enrichr/","GO_Molecular_Function_2021_", i, ".csv"))
#     write.csv(enriched[[2]], paste0(fig_output_path,"enrichr/","GO_Cellular_Component_2021_", i, ".csv"))
#     write.csv(enriched[[3]], paste0(fig_output_path,"enrichr/","GO_Molecular_Function_2021_", i, ".csv"))

#     # Make barplot of adjusted pvalues
#     BP <- enriched[[1]] %>%
#         arrange(Adjusted.P.value) %>%
#         head(10) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(Term, -Adjusted.P.value)))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," GO_Biological_Process_2021"))

#     CC <- enriched[[2]] %>%
#         arrange(Adjusted.P.value) %>%
#         head(10) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(Term, -Adjusted.P.value)))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," GO_Cellular_Component_2021"))

#     MF <- enriched[[3]] %>%
#         arrange(Adjusted.P.value) %>%
#         head(10) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(Term, -Adjusted.P.value)))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," GO_Molecular_Function_2021"))

#     kegg <- enriched[[4]] %>%
#         arrange(Adjusted.P.value) %>%
#         head(10) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(Term, -Adjusted.P.value)))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," KEGG_2021_Human"))

#     # Add to list so plots can be combined
#     ls_BP[[i]] <- BP
#     ls_CC[[i]] <- CC    
#     ls_MF[[i]] <- MF
#     ls_kegg[[i]] <- kegg
# }

# plot_grid(ls_BP[[1]],ls_BP[[2]],ls_BP[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","GO_Biological_Process_2021_all.png"), width=30, height=30, unit="cm")

# plot_grid(ls_CC[[1]],ls_CC[[2]],ls_CC[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","GO_Cellular_Component_2021_all.png"), width=30, height=30, unit="cm")

# plot_grid(ls_MF[[1]],ls_MF[[2]],ls_MF[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","GO_Molecular_Function_2021_all.png"), width=30, height=30, unit="cm")

# plot_grid(ls_kegg[[1]],ls_kegg[[2]],ls_kegg[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","KEGG_2021_Human_all.png"), width=30, height=30, unit="cm")

# plot_grid(ls_BP[[1]],ls_BP[[2]],ls_BP[[3]],ls_CC[[1]],ls_CC[[2]],ls_CC[[3]],ls_MF[[1]],ls_MF[[2]],ls_MF[[3]], 
#         nrow=3,ncol=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","all.png"), width=120, height=40, unit="cm", dpi =400)

# ###########################
# # fig 3 splice vs gene
# ###########################
if (!dir.exists(paste0(fig_output_path,"PS/withinType/"))){
dir.create(paste0(fig_output_path,"PS/withinType/"),
recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(paste0(fig_output_path,"PS/group_label/"))){
dir.create(paste0(fig_output_path,"PS/group_label/"),
recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(paste0(fig_output_path,"PS/main_label/"))){
dir.create(paste0(fig_output_path,"PS/main_label/"),
recursive = TRUE, showWarnings = TRUE)}

print(head(df_metadata))

# df_splice_ref_z %>%
#     arrange(desc(ratio)) %>%
#     head(10)

# df_splice_ref_z$high_ratio_2 <- NA

df_splice_ref_z_plot <- df_splice_ref_z %>%
    mutate(ratio = abs(splice_z/gene_z)) %>%
    mutate(Color = ifelse(ratio > 5, "ratio > 5", "ratio < 5")) %>%
    arrange(desc(ratio)) %>%
    filter(ratio != "NA")



# ############################
# # fig 3 zscore scatter plot
# ###############################
# Splice vs Gene________________________________________________________________
p_vs_gene <-  ggplot(aes(x=splice_z, y=gene_z, color = Color), data=df_splice_ref_z_plot) + 
    scale_color_manual(values=c("ratio > 5" = "red", "ratio < 5"="black")) +
    geom_point(size=.85, alpha = .6) +
    labs(colour = NULL, title= "", y = "Gene Expression Z-score", x = "Percent Spliced Z-score") +
    geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    # geom_text(
    #           label= df_splice_ref_z[["high_ratio_2"]],
    #           nudge_x = 0.01, nudge_y = 0.01,
    #           check_overlap =F, col = "black", size = 3
    #         ) +
    theme_classic() +
    theme(legend.position="right", 
            # legend.title=element_blank(),
        #   legend.title = element_text(size = 6), 
          legend.text = element_text(size = 8),
          axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12))  
        #   +
    # guides(color = guide_legend(nrow = 3)) 
    # +
    # ylim(-1, 1) 

ggsave(plot = p_vs_gene, dpi = 400,
    filename = paste0(fig_output_path, "/zscore_splice_ref_vs_gene.png"),  width=20, height= 8, unit="cm")


# ############################
# # fig 3 splice vs gene dist
# ################################

# plot_PS_gene <- function(junc,gene, comp_type, cell_type, df_PSorIR, tag ){
#     samples_main <- samples_other <- NA
#     print("____________________________________________")
#     print(cell_type)
#     # replace underscores, not sure why needed twice
#     cell_type_ <- sub("_", " ", cell_type)
#     cell_type_ <- sub("_", " ", cell_type_)

#     print(cell_type_)
#     print(comp_type)
#     print(junc)

#     if (comp_type == "main_label"){
#         samples_main <- df_metadata %>%
#                     filter(main_label == cell_type_) %>%
#                     pull(Run)
#         samples_other <- df_metadata %>%
#                     filter(main_label != cell_type_) %>%
#                     pull(Run)

#     } else if (comp_type == "group_label"){
#         samples_main <- df_metadata %>%
#                     filter(group_label == cell_type_) %>%
#                     pull(Run)
#         samples_other <- df_metadata %>%
#                     filter(group_label != cell_type_) %>%
#                     pull(Run)
#     } else {
#         # Within group comparison 
#         # get group label for the given label  
#         group <- df_metadata %>%
#                 filter(main_label == cell_type_) %>%
#                 pull(group_label)
#         group <- unique(group)
#         samples_main <- df_metadata %>%
#                      filter(main_label == cell_type_) %>%
#                      pull(Run)
#         samples_other <- df_metadata %>%
#                     filter(main_label != cell_type_) %>%
#                     filter(group_label == eval(group)) %>%
#                     pull(Run)
#     }

    
#     print(length(samples_main))
#     print(length(samples_other))

#     ##############
#     # PS Plot
#     #############

#     df_PS_junc <- df_PSorIR %>% 
#         filter(row.names(.) %in% c(junc)) %>%
#         t() %>%
#         as.data.frame() 

#     df_PS_junc_main <- df_PS_junc %>%
#         filter(row.names(.) %in% samples_main) %>%
#         mutate(type = cell_type)
#     df_PS_junc_other <- df_PS_junc %>%
#         filter(row.names(.) %in% samples_other) %>%
#         mutate(type = "other")

#     df_PS_junc_p <- rbind(df_PS_junc_main,df_PS_junc_other)
#     names(df_PS_junc_p) <- c("junc","type")

#     p1 <- ggplot(df_PS_junc_p, aes(x=type, y=junc)) +
#         geom_boxplot(outlier.shape = NA) +
#         theme_classic() +
#         geom_point(alpha = .5,position = position_jitter(seed = 1, width = 0.2)) +
#         # theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
#         theme(axis.title.x = element_blank()) + 
#         labs(y= tag, title = junc) 

#     ###################
#     # Expression Plot
#     ##################
#     df_gene <- df_exp %>% 
#         filter(row.names(.) %in% c(gene)) %>%
#         t() %>%
#         as.data.frame() 

#     df_gene_main <- df_gene %>%
#         filter(row.names(.) %in% samples_main) %>%
#         mutate(type = cell_type)
#     df_gene_other <- df_gene %>%
#         filter(row.names(.) %in% samples_other) %>%
#         mutate(type = "other")

#     df_gene_p <- rbind(df_gene_main,df_gene_other)
#     names(df_gene_p) <- c("junc","type")
#     p2 <- ggplot(df_gene_p,aes(x=type, y=junc)) +
#         geom_boxplot(outlier.shape = NA) +
#         theme_classic() +
#         geom_point(alpha = .5, position = position_jitter(seed = 1, width = 0.2)) +
#         # theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#         theme(axis.title.x = element_blank()) + 
#         labs(y= "log2(TPM+1)", title =gene) +
#         scale_y_continuous(
#             expand = c(0,0),
#             limits = c(0,  max(df_gene_p$junc) )) +
#         coord_cartesian(clip = "off") 


#     plot_grid(p1, p2)
#     ggsave(paste0(fig_output_path,tag,"/",comp_type,"/",junc,"_",gene,"_",cell_type,"_",comp_type,"_genePS.png"), width=15, height=8, unit="cm")

# }

# # for (row in 1:nrow(df_splice_ref_z)) {
# for (row in 1:100) {
#     junc <- as.character(df_splice_ref_z[row, "event"])
#     gene <- as.character(df_splice_ref_z[row, "overlapping"])
#     comp_type <- as.character(df_splice_ref_z[row, "group"])
#     cell_type <- as.character(df_splice_ref_z[row, "cell_type"])

#     plot_PS_gene(junc,gene, comp_type, cell_type , df_all_PS)

# }


#########################
# fig 4 IR vs gene
#########################
# if (!dir.exists(paste0(fig_output_path,"IR/withinType/"))){
# dir.create(paste0(fig_output_path,"IR/withinType/"),
# recursive = TRUE, showWarnings = TRUE)}
# if (!dir.exists(paste0(fig_output_path,"IR/group_label/"))){
# dir.create(paste0(fig_output_path,"IR/group_label/"),
# recursive = TRUE, showWarnings = TRUE)}
# if (!dir.exists(paste0(fig_output_path,"IR/main_label/"))){
# dir.create(paste0(fig_output_path,"IR/main_label/"),
# recursive = TRUE, showWarnings = TRUE)}

# print(head(df_metadata))
# print("----------")
# df_IR_ref_z  %>%
#     arrange(desc(ratio)) %>%
#     head(100)

# print(nrow(df_IR_ref_z))
# print(sum(is.na(df_IR_ref_z$gene_z)))
# print(sum(is.na(df_IR_ref_z$IR_z)))

# print(nrow(df_splice_ref_z))
# print(sum(is.na(df_splice_ref_z$gene_z)))
# print(sum(is.na(df_splice_ref_z$splice_z)))

# df_IR_ref_z_plot <- df_IR_ref_z %>%
#     mutate(ratio = abs(IR_z/gene_z)) %>%
#     mutate(Color = ifelse(ratio > 5, "ratio > 5", "ratio < 5")) %>%
#     arrange(desc(ratio)) %>%
#     filter(ratio != "NA")

# ############################
# # fig 4 IR zscore scatter plot
# ###############################
# # IR vs Gene________________________________________________________________
# IR_vs_gene <-  ggplot(aes(x=IR_z, y=gene_z, color = Color), data=df_IR_ref_z_plot) + 
#     scale_color_manual(values=c("ratio > 5" = "red", "ratio < 5"="black")) +
#     geom_point(size=.85, alpha = .6) +
#     labs(colour = NULL, title= "", y = "Gene Expression Z-score", x = "Intron Retention Z-score") +
#     geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
#     geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
#     theme_classic() +
#     theme(legend.position="right", 
#             # legend.title=element_blank(),
#         #   legend.title = element_text(size = 6), 
#           legend.text = element_text(size = 8),
#           axis.title.x=element_text(size=12),
#           axis.title.y=element_text(size=12))  
#         #   +
#     # guides(color = guide_legend(nrow = 3)) 
#     # +
#     # ylim(-1, 1) 

# ggsave(plot = IR_vs_gene, dpi = 400,
#     filename = paste0(fig_output_path, "/zscore_IR_ref_vs_gene.png"),  width=20, height= 8, unit="cm")

# for (row in 1:100) {
#     junc <- as.character(df_IR_ref_z[row, "event"])
#     gene <- as.character(df_IR_ref_z[row, "overlapping"])
#     comp_type <- as.character(df_IR_ref_z[row, "group"])
#     cell_type <- as.character(df_IR_ref_z[row, "cell_type"])

#     plot_PS_gene(junc,gene, comp_type, cell_type , df_IR_table, "IR")

# }

# #############################
# # comparing sets 
# ################################

# # Get union of sets in both to iterate through 
# ls_all_sets <- union(df_splice_set$set, df_gene_set$set)
# # print(length(ls_all_sets))

# # Get intersection of spliced genes and genes
# gene_intersection <- intersect(unique(df_splice_set$overlapping), unique(df_gene_set$X))
# # print(gene_intersection)

# cat("\n Number of common genes in gene and splice sets: \n")
# print(length(gene_intersection))

# cat("\n Number of unique genes in splice set: \n")
# print(length(unique(df_splice_set$overlapping)))

# cat("\n Number of unique junctions in splice set: \n")
# print(length(unique(df_splice_set$event)))

# cat("\n Number of unique genes in gene set: \n")
# print(length(unique(df_gene_set$X)))
# cat("\n")

# ###############################
# # Compare on set level 
# ###############################
# df_counts <- data.frame(matrix(ncol = 3, nrow = length(ls_all_sets)))
# rownames(df_counts) <- ls_all_sets
# colnames(df_counts) <- list("gene", "splice","IR")

# for (s in ls_all_sets){
#   print(s)

#   # Counting 
#   df_splice <- df_splice_set %>% 
#     filter(set == s) 
#   df_IR <- df_IR_set %>% 
#     filter(set == s) 
#   df_gene <- df_gene_set %>% 
#     filter(set == s) 
#   df_counts[s, "splice"] <- nrow(df_splice)
#   df_counts[s, "gene"] <- nrow(df_gene)
#   df_counts[s, "IR"] <- nrow(df_IR)

#   # Comparing 
# #   gene_intersection <- intersect(df_splice$overlapping, df_gene$X)
# #   print(length(gene_intersection))
# #   df_counts[s, "gene and splice"] <- length(gene_intersection)
# }

# print(head(df_counts))

# df_counts %>% 
#   pivot_longer(names_to = "variable", values_to = "value", cols = everything()) %>%
#   ggplot(aes(x = fct_inorder(variable), y = value, color = variable)) +
#     geom_jitter(alpha = 0.5, size = 2) +
#     theme_classic() +
#     theme(legend.position = "none") +
#     labs(x = "", y = "marker count per set") +
#     scale_color_manual(values=c("gene"="orange","splice"="skyblue","IR"="darkgreen"))

# ggsave(paste0(fig_output_path, "count_per_set_plot.png"), width=6, height=10, dpi=300, units=c("cm"))

# # df_counts %>% 
# #   tibble::rownames_to_column("set") %>%
# #   arrange(desc(splice))

# df_counts %>% 
#   tibble::rownames_to_column("set") %>%
#   filter(gene < 100) %>%
#   arrange(desc(gene))

