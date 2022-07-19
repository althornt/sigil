#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(uwot)
library(purrr)
library(ComplexHeatmap)
library(UpSetR)
library(cowplot)
library(enrichR)
library(GenomicRanges)
library(valr)
library(rGREAT)
library(ggpubr)
library(ggrepel)

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

# group <- ggplot(df_metadata, aes(x=group_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#           axis.title.y = element_text(size=8),
#           legend.position="None", 
#             axis.title.x = element_blank()) +
#   scale_fill_manual(values=c("#999999", "black")) +
#     ylab("Number of samples")

# ggsave(paste0(fig_output_path,"metadata_bar_group.png"), width=10, height=5.5, unit="cm")

# ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#     legend.position="None",
#     axis.title.x = element_blank()) +
#     scale_fill_manual(values=c("#999999", "black"))+
#     scale_y_continuous(breaks= scales::pretty_breaks())

# ggsave(paste0(fig_output_path,"metadata_bar_main.png"), width=10, height=15, unit="cm")

# main<-ggplot(df_metadata, aes(x=main_label, fill=data_source)) +
#   geom_bar(position="stack") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.title.y = element_text(size=8),
#             legend.position="top",   
#             axis.title.x = element_blank(),
#             legend.title=element_text(size=8),
#             legend.text=element_text(size=8))+
#             # legend.box="vertical", 
#             # legend.margin=margin())+
#     # guides(color = guide_legend(override.aes = list(size=1)))+
#     scale_fill_manual("data source", values=c("#999999", "black"))+
#     scale_y_continuous(breaks= scales::pretty_breaks()) +
#     ylab("Number of samples")

# ggsave(paste0(fig_output_path,"metadata_bar_main_legend.png"), width=10, height=15, unit="cm")

# # ggpubr::ggarrange(main, group, # list of plots
# #                 #   labels = "AUTO", # labels
# #                   common.legend = T, # COMMON LEGEND
# #                   legend = "top", # legend position
# #                   align = "hv", # Align them both, horizontal and vertical
# #                   nrow = 2)  # number of rows

# pg <- plot_grid(main, group, nrow=2, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"metadata_bar_combined.png"), width=10, height=20, unit="cm")
# #############################
# # Fig 2A comparing sets 
# ################################

# Get union of sets in both to iterate through 
ls_all_sets <- union(df_splice_set$set, df_gene_set$set)
# print(length(ls_all_sets))

# Get intersection of spliced genes and genes
gene_intersection <- intersect(unique(df_splice_set$overlapping), unique(df_gene_set$X))
# print(gene_intersection)

cat("\n Number of common genes in gene and splice sets: \n")
print(length(gene_intersection))

cat("\n Number of unique genes in splice set: \n")
print(length(unique(df_splice_set$overlapping)))

cat("\n Number of unique junctions in splice set: \n")
print(length(unique(df_splice_set$event)))

cat("\n Number of unique genes in gene set: \n")
print(length(unique(df_gene_set$X)))
cat("\n")

# ###############################
# # Fig 2A Compare on set level 
# ###############################
df_counts <- data.frame(matrix(ncol = 3, nrow = length(ls_all_sets)))
rownames(df_counts) <- ls_all_sets
colnames(df_counts) <- list("Gene", "Splice","IR")

for (s in ls_all_sets){
  print(s)

  # Counting 
  df_splice <- df_splice_set %>% 
    filter(set == s) 
  df_IR <- df_IR_set %>% 
    filter(set == s) 
  df_gene <- df_gene_set %>% 
    filter(set == s) 
  df_counts[s, "Splice"] <- nrow(df_splice)
  df_counts[s, "Gene"] <- nrow(df_gene)
  df_counts[s, "IR"] <- nrow(df_IR)

  # Comparing 
#   gene_intersection <- intersect(df_splice$overlapping, df_gene$X)
#   print(length(gene_intersection))
#   df_counts[s, "gene and splice"] <- length(gene_intersection)
}

print(head(df_counts))

#     cell_type_ <- sub("_", " ", cell_type)


df_counts %>% 
    rownames_to_column("set") %>%
    rowwise() %>%
  pivot_longer(!set, names_to = "variable", values_to = "value") %>%
  mutate(setlabel = "") %>%
  mutate(setlabel = ifelse(((value < 100)&(variable=="Gene")), set, setlabel)) %>%
  # mutate(setlabel = sub("_", " ", setlabel)) %>%
  # mutate(setlabel = sub("_", " ", setlabel)) %>% #needed 3
  ggplot(aes(x = fct_inorder(variable), y = value, color = variable)) +
    geom_jitter(alpha = 0.5, size = 1.5) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(size=10), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10)) +
    labs(x = "", y = "Marker count per set") +
    scale_color_manual(values=c("Gene"="orange","Splice"="skyblue","IR"="darkgreen"))
    # +
#     coord_cartesian(clip = "off")+
#     geom_text_repel(aes(label = setlabel), size = 3.5,
#         xlim = c(NA, Inf), ylim = c(NA, Inf),
#         # min.segment.length = 0, 
#       force        = 2,
#     # # nudge_x      = 2,
#     # direction    = "x",
#     hjust        = 1, segment.square = FALSE, 
#     # segment.size = 0.2,
#     # box.padding = 0.3
#   )+
#   scale_x_discrete(
#   # breaks = 1:2, labels = c("Gene", "Splice", "Intron Retention"),
#   expand = expansion(mult = 1)
# )
#     # # Repel away from the left edge, not from the right.
#     # xlim = c(NA, Inf),
#     # # Do not repel from top or bottom edges.
#     # ylim = c(-Inf, Inf), fill = "white") +

ggsave(paste0(fig_output_path, "count_per_set_plot.png"), width=6, height=8, dpi=300, units=c("cm"))

# # df_counts %>% 
# #   tibble::rownames_to_column("set") %>%
# #   arrange(desc(splice))

# df_counts %>% 
#   tibble::rownames_to_column("set") %>%
#   filter(Gene < 100) %>%
#   arrange(desc(Gene))

# #######################################################
# # Fig2B UpSet plot comparing type of sets gene usage
# ######################################################
# Genes affected by splice, IR, and gene
ls_gene_ref_genes <- unique(df_gene_set$X)
ls_splice_ref_genes <- unique(df_splice_set$overlapping)
ls_IR_ref_genes <- unique(df_IR_set$overlapping)

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

# df_splice_set %>%
#     filter(overlapping %in% ls_ref_genes_splice_IR) %>%
#     select(event, overlapping, cell_type, group)


# df_IR_set %>%
#     filter(overlapping %in% ls_ref_genes_splice_IR) %>%
#     select(event, overlapping, cell_type, group)

# for (i in ls_IR_ref_genes){
#     cat("\n",i)
# }


print("Number of gene ref genes:")
print(length(ls_gene_ref_genes))
print("Number of splice ref genes:")
print(length(ls_splice_ref_genes))
print("Number of IR ref genes:")
print(length(ls_IR_ref_genes))

# Gene level 
listInput <- list("Gene" = ls_gene_ref_genes,
                 "Splice" = ls_splice_ref_genes, 
                 "Intron Retention" = ls_IR_ref_genes)

# plotObject <- upset(fromList(listInput), order.by = "freq", text.scale = c(1.6, 1.6, 1.6, 1.6, 2, 1.3),sets.bar.color=c("orange","skyblue","darkgreen"))
# png(file= paste0(fig_output_path, "/upsetplot_ref_genes_vs_splice_vs_IR.png"),
#     width = 20,
#     height    = 14,
#     units     = "cm",
#     res       = 1200)
# print(plotObject)
# dev.off()

# # Junction level 
# listInput <- list(
#                  "Splice" = unique(df_splice_set$event), 
#                  "Intron Retention" = unique(df_IR_set$event))

# plotObject <- upset(fromList(listInput), order.by = "freq", text.scale = c(1.6, 1.6, 1.6, 1.6, 2, 1.3),sets.bar.color=c("skyblue","darkgreen"))
# png(file= paste0(fig_output_path, "/upsetplot_junctions_splice_vs_IR.png"),
#     width = 20,
#     height    = 14,
#     units     = "cm",
#     res       = 1200)
# print(plotObject)
# dev.off()

# ##########################################################
# # Fig 3 compare junctions to UCSC Alt track alt promoter
# ############################################################
# # Add bed cols to splice ref 
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

# # Get number to randomly sample from PS file 
# num_unique_ref_juncs <- length(unique(df_splice_set$event))
# print(length(unique(num_unique_ref_juncs)))

# # Randonmly sample from the file with all junctions 
# ls_random_sample <- sample(rownames(df_all_PS), num_unique_ref_juncs, replace=F)

# print(head(ls_random_sample))
# print(length(ls_random_sample))
# print(length(unique(ls_random_sample)))

# print("#############")

# # Format random junctions into bed 
# df_random_juncs <- as.data.frame(ls_random_sample) %>%
#     dplyr::rename(event = ls_random_sample )  %>%
#     rowwise() %>%
#     mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1]))  %>%
#     mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
#     mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
#     mutate(strand = strsplit(event, ":")[[1]][3]) %>%
#     as.data.frame() 


# print(head(df_random_juncs))
# print(dim(df_random_juncs))


# # Format all junctions into bed 
# ls_all_juncs <- rownames(df_all_PS)
# print(length(ls_all_juncs))
# ls_all_juncs <- ls_all_juncs[!(ls_all_juncs %in% df_splice_set$event)]
# print(length(ls_all_juncs))

# df_all_background_juncs <- as.data.frame(ls_all_juncs) %>%
#     dplyr::rename(event = ls_all_juncs )  %>%
#     rowwise() %>%
#     mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1]))  %>%
#     mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
#     mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
#     mutate(strand = strsplit(event, ":")[[1]][3]) %>%
#     as.data.frame() 


# print(head(df_all_background_juncs))
# print(dim(df_all_background_juncs))


# # Set up alt promoter df
# df_alt_events_promoter <- df_alt_events %>%
#     filter(name=="altPromoter") %>%
#     select("chrom","start","end","name","strand")


# compAltProm <- function(df_bed ){

#     # Find intersection of splice ref junctions with alt promoters; 
#     # if they intersect/overlap it suggests the alt promoter its intersected with 
#     # is being skipped
#     df_alt_promoter_splice_ref_int <- bed_intersect(df_bed, df_alt_events_promoter, 
#             suffix = c(".splice_ref", ".alt_pro"))  %>%
#             as.data.frame() 
#     print(head(df_alt_promoter_splice_ref_int))

#     # recurring_int <- df_alt_promoter_splice_ref_int %>%
#     #     group_by(event.splice_ref) %>%
#     #     count() %>%
#     #     arrange(desc(n)) %>%
#     #     as.data.frame()

#     # print(head(recurring_int, n=40))
#     # print(dim(df_alt_promoter_splice_ref_int))
#     # print(dim(df_bed))

#     # Find which junctions within 200 of an alt promoter 
#     # but do not overlap with an alt promoter, since thats accounted for in the other value
#     df_alt_promoter_splice_ref_closest <- bed_closest(df_bed, df_alt_events_promoter, 
#             suffix = c(".splice_ref", ".alt_pro"))  %>%
#             as.data.frame() %>%
#             filter(abs(.dist) <= 200) %>%
#             filter(.overlap == 0) %>%
#             as.data.frame()  
#     head(df_alt_promoter_splice_ref_closest, n=20)

#     # Count results | need unique because the input df has repeat events for every set its present in
#     unique_splice_juncs <- unique(df_bed$event)
#     unique_splice_juncs_alt_pro_int <- unique(df_alt_promoter_splice_ref_int$event.splice_ref)
#     unique_splice_juncs_alt_pro_close <- unique(df_alt_promoter_splice_ref_closest$event.splice_ref)
#     unique_splice_juncs_alt_pro_either <- union(df_alt_promoter_splice_ref_int$event.splice_ref, 
#                                                 df_alt_promoter_splice_ref_closest$event.splice_ref )

#     cat("\n Percent of splice ref junctions that contain an alt promoter: \n")
#     perc_contain_alt <- length(unique_splice_juncs_alt_pro_int)/length(unique_splice_juncs)*100 
#     cat(perc_contain_alt)
#     cat("\n")

#     cat("\n Percent of splice ref junctions that are 200 bp downstream of an alt promoter: \n")
#     perc_downstream_alt <- length(unique_splice_juncs_alt_pro_close)/length(unique_splice_juncs)*100
#     cat(perc_downstream_alt )
#     cat("\n")

#     cat("\n Percent of splice ref junctions that are either: \n")
#     percent_alt <-length(unique_splice_juncs_alt_pro_either)/length(unique_splice_juncs)*100 
#     cat(percent_alt)
#     cat("\n")
#     # Make alt pro calc res into df
#     df_res_perc <- t(data.frame(list("Intersects"=perc_contain_alt,
#                                     "Downstream"=perc_downstream_alt, 
#                                     "Alternative Promoter"=percent_alt )))
#     print(df_res_perc)
#     df_res_counts <- t(data.frame(list("Intersects"= length(unique_splice_juncs_alt_pro_int),
#                                     "Downstream"=length(unique_splice_juncs_alt_pro_close),
#                                      "Alternative Promoter"=length(unique_splice_juncs_alt_pro_either) )))
#     print(df_res_counts)
#     # Add alt promoter info the input df 
#     df_bed <- df_bed %>%
#         mutate(contains_alt_pro = ifelse((event %in% unique_splice_juncs_alt_pro_int), TRUE, FALSE)) %>%
#         mutate(downstream_alt_pro = ifelse((event %in% unique_splice_juncs_alt_pro_close), 1, 0)) %>%
#         mutate(contains_and_downstream_alt_pro = ifelse(((event %in% unique_splice_juncs_alt_pro_int) & (event %in% unique_splice_juncs_alt_pro_close)), 1, 0)) %>%
#         mutate(AltPromoter = ifelse(((event %in% unique_splice_juncs_alt_pro_int) | (event %in% unique_splice_juncs_alt_pro_close)), TRUE, FALSE)) 

#     return(list(df_bed, df_res_perc, df_res_counts))

# }


# ls_alt_res_splice_sets <- compAltProm(df_splice_set)
# df_splice_set <- ls_alt_res_splice_sets[[1]] 
# head(df_splice_set)

# # print("~~~~~~~~~~~")
# ls_alt_res_random <- compAltProm(df_random_juncs)
# df_random_juncs <- ls_alt_res_random[[1]]
# head(df_random_juncs)

# print("~~~~~~~~~~~")
# ls_alt_res_background <- compAltProm(df_all_background_juncs)
# df_all_background <- ls_alt_res_background[[1]]
# head(df_all_background)

# # Make plot comparing alt promoter count in splice sets to random junctions
# df_alt_pro_res <- cbind(ls_alt_res_splice_sets[[2]], ls_alt_res_random[[2]], ls_alt_res_background[[2]])
# colnames(df_alt_pro_res) <- list("splice_ref", "random_sample", "background")


# # Fishers exact test 

# ##                      Non-AFE      AFE
# ## Splice ref           4.5           4.5
# ## Non splice ref       2.5             2.5
# df_alt_pro_res_count <- cbind(ls_alt_res_splice_sets[[3]], ls_alt_res_background[[3]])
# colnames(df_alt_pro_res_count) <- list("splice_ref", "background")
# print(df_alt_pro_res_count)


# # Building contingency table 
# AltPro_splice_ref <- df_alt_pro_res_count["Alternative.Promoter","splice_ref"]
# AltPro_splice_background <- df_alt_pro_res_count["Alternative.Promoter","background"]

# NotAltPro_splice_ref <- num_unique_ref_juncs - AltPro_splice_ref
# NotAltPro_splice_background <- length(ls_all_juncs) - AltPro_splice_background

# print(paste0(AltPro_splice_ref, AltPro_splice_background, NotAltPro_splice_ref, NotAltPro_splice_background, sep= " "))

# df_fishers <- data.frame("nonAFE" = c(NotAltPro_splice_ref, NotAltPro_splice_background), 
#                 "AFE" = c(AltPro_splice_ref, AltPro_splice_background), 
#                 row.names = c("spliceRef", "background"))


# # Run Fishers test
# print(df_fishers)

# res_fishers <- stats::fisher.test(df_fishers)
# print(res_fishers)
# print(res_fishers$p.value)

# df_alt_pro_res_perc <- cbind(ls_alt_res_splice_sets[[2]], ls_alt_res_background[[2]])
# colnames(df_alt_pro_res_perc) <- list("Splice set junctions", "All junctions")
# print(df_alt_pro_res_perc)


# df_alt_pro_res_perc <- df_alt_pro_res_perc %>% 
#     as.data.frame() %>%
#     rownames_to_column( var= "type")
# print(df_alt_pro_res_perc)

# df_alt_pro_res_perc <- pivot_longer(df_alt_pro_res_perc,  values_to = "Percent",  
#                                   names_to = "name",  
#                                   cols = c("Splice set junctions", "All junctions")) 
# print(df_alt_pro_res_perc)
# df_alt_pro_res_perc$type <- sub("[.]", " ", df_alt_pro_res_perc$type)
# df_alt_pro_res_perc$type <- sub("[.]", " ", df_alt_pro_res_perc$type)

# # Number of cars in each class:
# df_alt_pro_res_perc %>% 
#     ggplot(aes(x=type, y=Percent, fill=name)) +
#     geom_bar(stat = 'identity', position = "dodge") +
#     theme_classic() +
#     scale_fill_manual("name", values=c("#999999", "black"))+
#     theme(legend.title= element_blank(), legend.text = element_text(size=12),legend.position = "bottom",
#      axis.title.x=element_blank(), 
#         axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
#     labs(title = paste0("Fishers test p-value = ", signif(res_fishers$p.value, 4))) 


# ggsave(paste0(fig_output_path,"alt_promoter.png"), width=20, height=10, unit="cm", dpi = 300)




# ##########################################################
# # Fig 3 rGREAT splice ref and IR ref GO terms
# ############################################################

# # Run Great on Splice set junctions 
# job_PS = submitGreatJob(makeGRangesFromDataFrame(df_splice_set),
#                         species= "hg38", version = "4")

# # Get tables from Run
# tb_PS = getEnrichmentTables(job_PS, category = c("GO"))

# plotGREAT <- function(tb, GO_type, ggtitle){
#     # Make barplot of adjusted pvalues for top terms
#     plt <- tb[[GO_type]] %>%
#       filter((Binom_Adjp_BH <= 0.05) &(Hyper_Adjp_BH <= 0.05 )) %>%
#       filter((Binom_Fold_Enrichment >= 2) &(Hyper_Fold_Enrichment >= 2 )) %>%
#       arrange(Binom_Adjp_BH) %>%
#       mutate(Binom_Adjp_BH_log10 = -log10(Binom_Adjp_BH)) %>%
#         head(15) %>%
#         ggplot(aes(x = Binom_Adjp_BH_log10, y = reorder(name, Binom_Adjp_BH_log10), fill =Binom_Fold_Enrichment))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 22),
#                   axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20),
#                   plot.title = element_text(size=25),
#                   legend.title = element_text(size = 20),legend.text = element_text(size = 20),
#                   legend.position = c(0.8, 0.1)) +
#            labs(title = paste0(ggtitle), x = "-log10(Adjusted P-value)", fill = "Fold Enrichment")+
#           scale_y_discrete(labels = function(x) str_wrap(x, width =40))

#     ggsave(paste0(fig_output_path,"rGREAT_", paste(ggtitle,sep="_"),".png"), 
#                 width=30, height=25, unit="cm", dpi = 300)

#   return(plt)
# }

# # Make plots for PS great results 
# great_PS_BP <- plotGREAT(tb_PS, "GO Biological Process", "Splice GO Biological Process" )
# great_PS_MF <-plotGREAT(tb_PS, "GO Molecular Function", "Splice GO Molecular Function" )

# # Combined PS plot 
# plot_grid(great_PS_BP, great_PS_MF, nrow=1, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"rGREAT_GO_MF_BP_all_splice_ref.png"), 
#             width=70, height=35, unit="cm", dpi=300)

# # Format IR table to add the needed bed columns
# df_IR_set <- df_IR_set %>%
#   rowwise() %>%
#     mutate(chrom = paste0("chr",strsplit(event, ":")[[1]][1])) %>%
#     mutate(start = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][1] )) %>%
#     mutate(end = as.numeric(strsplit(strsplit(event, ":")[[1]][2], "-")[[1]][2] ))    %>%
#     mutate(strand = strsplit(event, ":")[[1]][3]) %>%
#     as.data.frame() 

# # Run Great on IR sigil junctions
# job_IR = submitGreatJob(makeGRangesFromDataFrame(df_IR_set),
#                         species= "hg38", version = "4")

# # Get tables from Run
# tb_IR = getEnrichmentTables(job_IR, category = c("GO"))

# # Make plots for IR great results 
# great_IR_BP <- plotGREAT(tb_IR, "GO Biological Process", "Intron Retention GO Biological Process" )
# great_IR_MF <-plotGREAT(tb_IR, "GO Molecular Function", "Intron Retention GO Molecular Function" )

# plot_grid(great_IR_BP, great_IR_MF, nrow=1, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"rGREAT_GO_MF_BP_all_IR_ref.png"), 
#             width=70, height=30, unit="cm", dpi=300)

# # Combined PS and IR GREAT plots
# plot_grid(great_PS_BP, great_PS_MF, great_IR_BP, great_IR_MF,
#              nrow=2, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"rGREAT_GO_MF_BP_all_PS_IR_ref.png"), width=75, height=60, unit="cm")

# # ##########################################################
# # # Fig 2 Enrichr on all genes from each set type
# # ############################################################
# if (!dir.exists(paste0(fig_output_path,"enrichr"))){
# dir.create(paste0(fig_output_path,"enrichr"),
# recursive = TRUE, showWarnings = TRUE)}


# dbs <- c("GO_Biological_Process_2021",
#         "GO_Cellular_Component_2021", 
#         "GO_Molecular_Function_2021",
#         "KEGG_2021_Human")

# ls_sigil_types <- list("Gene"=ls_gene_ref_genes,
#                     "Splice"=ls_splice_ref_genes,
#                     "Intron Retention"=ls_IR_ref_genes)
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
#       rowwise() %>%
#         mutate(term_clean = strsplit(Term, "[(]" )[[1]][1]) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(term_clean, -Adjusted.P.value)))  + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," GO BP"))

#     CC <- enriched[[2]] %>%
#         arrange(Adjusted.P.value) %>%
#         head(10) %>%
#       rowwise() %>%
#         mutate(term_clean = strsplit(Term, "[(]" )[[1]][1]) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(term_clean, -Adjusted.P.value)))   + 
#                     geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," GO CC"))

#     MF <- enriched[[3]] %>%
#         arrange(Adjusted.P.value) %>%
#         head(10) %>%
#       rowwise() %>%
#         mutate(term_clean = strsplit(Term, "[(]" )[[1]][1]) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(term_clean, -Adjusted.P.value)))   + 
#                     geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," GO MF"))

#     kegg <- enriched[[4]] %>%
#         arrange(Adjusted.P.value) %>%
#               head(10) %>%
#       rowwise() %>%
#         mutate(term_clean = strsplit(Term, "[(]" )[[1]][1]) %>%
#         ggplot(aes(x = -log(Adjusted.P.value), y = reorder(term_clean, -Adjusted.P.value)))   + 
#             geom_bar(stat = 'identity') +
#             theme_classic() +
#             theme(axis.title.y=element_blank(), axis.text = element_text(size = 15), plot.title = element_text(size=15)) +
#             labs(title = paste0(i," KEGG"))

#     # Add to list so plots can be combined
#     ls_BP[[i]] <- BP
#     ls_CC[[i]] <- CC    
#     ls_MF[[i]] <- MF
#     ls_kegg[[i]] <- kegg
# }

# plot_grid(ls_BP[[1]],ls_BP[[2]],ls_BP[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","GO_Biological_Process_2021_all.png"), width=30, height=30, unit="cm", dpi =400)

# plot_grid(ls_CC[[1]],ls_CC[[2]],ls_CC[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","GO_Cellular_Component_2021_all.png"), width=30, height=30, unit="cm",dpi =400)

# plot_grid(ls_MF[[1]],ls_MF[[2]],ls_MF[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","GO_Molecular_Function_2021_all.png"), width=30, height=30, unit="cm",dpi =400)

# plot_grid(ls_kegg[[1]],ls_kegg[[2]],ls_kegg[[3]], nrow=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","KEGG_2021_Human_all.png"), width=30, height=30, unit="cm",dpi =400)

# plot_grid(ls_BP[[1]],ls_BP[[2]],ls_BP[[3]],ls_CC[[1]],ls_CC[[2]],ls_CC[[3]],ls_MF[[1]],ls_MF[[2]],ls_MF[[3]], 
#         nrow=3,ncol=3, align = 'v', axis = 'l')
# ggsave(paste0(fig_output_path,"enrichr/","all.png"), width=120, height=40, unit="cm", dpi =400)

# # ###########################
# # # fig 3 splice vs gene
# # ###########################
# if (!dir.exists(paste0(fig_output_path,"PS/withinType/"))){
# dir.create(paste0(fig_output_path,"PS/withinType/"),
# recursive = TRUE, showWarnings = TRUE)}
# if (!dir.exists(paste0(fig_output_path,"PS/group_label/"))){
# dir.create(paste0(fig_output_path,"PS/group_label/"),
# recursive = TRUE, showWarnings = TRUE)}
# if (!dir.exists(paste0(fig_output_path,"PS/main_label/"))){
# dir.create(paste0(fig_output_path,"PS/main_label/"),
# recursive = TRUE, showWarnings = TRUE)}

# print(head(df_metadata))

# # df_splice_ref_z %>%
# #     arrange(desc(ratio)) %>%
# #     head(10)

# df_splice_ref_z$high_ratio_2 <- NA

df_splice_ref_z_plot <- df_splice_ref_z %>%
    mutate(ratio = abs(splice_z/gene_z)) %>%
    mutate(Color = ifelse(ratio > 5, "ratio > 5", "ratio < 5")) %>%
    arrange(desc(ratio)) %>%
    filter(ratio != "NA")

# print(head(df_splice_set[c("event","contains_alt_pro","downstream_alt_pro","contains_and_downstream_alt_pro","AltPromoter")]))
# print(head(df_splice_ref_z_plot))

# # Merge alt promoter cols into splice df
# df_splice_ref_z_plot <- merge(df_splice_set[c("event","contains_alt_pro","downstream_alt_pro","contains_and_downstream_alt_pro","AltPromoter")],
#             df_splice_ref_z_plot,by="event")

# print(dim(df_splice_ref_z_plot))

# print("_________________________________")
# print(head(df_splice_ref_z_plot))
# print(dim(df_splice_ref_z_plot))

# # count how many are missing alt pro info 

# # ############################
# # # fig 3 zscore scatter plot
# # ###############################
# Splice vs Gene________________________________________________________________
p_vs_gene <-  ggplot(aes(x=splice_z, y=gene_z, color = Color), data=df_splice_ref_z_plot) + 
    scale_color_manual(values=c("ratio > 5" = "red", "ratio < 5"="black")) +
    geom_point(size= 1, alpha = .8) +
    labs(colour = NULL, title= "", y = "Gene Expression Z-score", x = "Percent Spliced Z-score") +
    geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    # geom_text(
    #           label= df_splice_ref_z[["high_ratio_2"]],
    #           nudge_x = 0.01, nudge_y = 0.01,
    #           check_overlap =F, col = "black", size = 3
    #         ) +
    theme_classic() +
    theme(
          legend.text = element_text(size = 10),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          legend.position = c(0.05, .9)
          )  +
    xlim(-6,6)
    # guides(color = guide_legend(nrow = 3)) 
    # +
    # ylim(-1, 1) 

ggsave(plot = p_vs_gene, dpi = 400,
    filename = paste0(fig_output_path, "/zscore_splice_ref_vs_gene.png"),  width=20, height= 8, unit="cm")


PS_z_hist <- ggplot(df_splice_ref_z_plot, aes(x=splice_z)) + 
    geom_histogram(bins = 100, colour="black", fill="black")+
    theme_classic()+
    xlim(-6, 6) +
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    labs(x = "Percent Spliced Z-score") +
    theme(   axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12)) 


ggsave(plot = PS_z_hist, dpi = 400,
    filename = paste0(fig_output_path, "/zscore_splice_hist.png"),  width=20, height= 8, unit="cm")

plot_grid(PS_z_hist, p_vs_gene, ncol = 1, rel_heights = c(.35, 1))
ggsave(file.path(fig_output_path,
                paste0("zscore_splice_ref_vs_gene_with_hist.png")),
        device = "png",
        width = 12, height = 5,
        dpi = 300, bg="white")
# Count number of positive and negative PS z-scores
table(sign(df_splice_ref_z_plot$splice_z))

# # coloring Alt promoter
# # Splice vs Gene________________________________________________________________
# p_vs_gene <-  ggplot(aes(x=splice_z, y=gene_z, color = AltPromoter), data=df_splice_ref_z_plot) + 
#     # scale_color_manual(values=c("ratio > 5" = "red", "ratio < 5"="black")) +
#     geom_point(size=.85, alpha = .6) +
#     labs(colour = NULL, title= "", y = "Gene Expression Z-score", x = "Percent Spliced Z-score") +
#     geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
#     geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
#     # geom_text(
#     #           label= df_splice_ref_z[["high_ratio_2"]],
#     #           nudge_x = 0.01, nudge_y = 0.01,
#     #           check_overlap =F, col = "black", size = 3
#     #         ) +
#     theme_classic() +
#     theme(legend.position="right", 
#           legend.title = element_text(size = 6), 
#           legend.text = element_text(size = 8),
#           axis.title.x=element_text(size=12),
#           axis.title.y=element_text(size=12))  

# ggsave(plot = p_vs_gene, dpi = 400,
#     filename = paste0(fig_output_path, "/zscore_splice_ref_vs_gene_AltProm.png"),  width=20, height= 8, unit="cm")


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
if (!dir.exists(paste0(fig_output_path,"IR/withinType/"))){
dir.create(paste0(fig_output_path,"IR/withinType/"),
recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(paste0(fig_output_path,"IR/group_label/"))){
dir.create(paste0(fig_output_path,"IR/group_label/"),
recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(paste0(fig_output_path,"IR/main_label/"))){
dir.create(paste0(fig_output_path,"IR/main_label/"),
recursive = TRUE, showWarnings = TRUE)}

print(head(df_metadata))
print("----------")
df_IR_ref_z  %>%
    arrange(desc(ratio)) %>%
    head(100)

print(nrow(df_IR_ref_z))
print(sum(is.na(df_IR_ref_z$gene_z)))
print(sum(is.na(df_IR_ref_z$IR_z)))

print(nrow(df_splice_ref_z))
print(sum(is.na(df_splice_ref_z$gene_z)))
print(sum(is.na(df_splice_ref_z$splice_z)))

df_IR_ref_z_plot <- df_IR_ref_z %>%
    mutate(ratio = abs(IR_z/gene_z)) %>%
    mutate(Color = ifelse(ratio > 5, "ratio > 5", "ratio < 5")) %>%
    arrange(desc(ratio)) %>%
    filter(ratio != "NA")

# ############################
# # fig 4 IR zscore scatter plot
# ###############################
# IR vs Gene________________________________________________________________
IR_vs_gene <-  ggplot(aes(x=IR_z, y=gene_z, color = Color), data=df_IR_ref_z_plot) + 
    scale_color_manual(values=c("ratio > 5" = "red", "ratio < 5"="black")) +
    geom_point(size=1, alpha = .8) +
    labs(colour = NULL, title= "", y = "Gene Expression Z-score", x = "Intron Retention Z-score") +
    geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    theme_classic() +
    theme(
          legend.text = element_text(size = 10),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          legend.position = c(0.05, .9)
          )  +
    guides(color = guide_legend(nrow = 2))+ 
    xlim(-6, 6) 

ggsave(plot = IR_vs_gene, dpi = 400,
    filename = paste0(fig_output_path, "/zscore_IR_ref_vs_gene.png"),  width=20, height= 8, unit="cm")



# for (row in 1:100) {
#     junc <- as.character(df_IR_ref_z[row, "event"])
#     gene <- as.character(df_IR_ref_z[row, "overlapping"])
#     comp_type <- as.character(df_IR_ref_z[row, "group"])
#     cell_type <- as.character(df_IR_ref_z[row, "cell_type"])

#     plot_PS_gene(junc,gene, comp_type, cell_type , df_IR_table, "IR")

# }

IR_PS_z_hist <- ggplot(df_IR_ref_z_plot, aes(x=IR_z)) + 
    geom_histogram(bins = 100, colour="black", fill="black")+
    theme_classic()+
    xlim(-6, 6) +
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    labs(x = "Intron Retention Z-score") +
    theme(   axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12)) 


ggsave(plot = IR_PS_z_hist, dpi = 400,
    filename = paste0(fig_output_path, "/zscore_IR_hist.png"),  width=20, height= 8, unit="cm")

plot_grid(IR_PS_z_hist, IR_vs_gene, ncol = 1, rel_heights = c(.35, 1))
ggsave(file.path(fig_output_path,
                paste0("zscore_IR_ref_vs_gene_with_hist.png")),
        device = "png",
        width = 12, height = 5,
        dpi = 300, bg="white")


  
# Count number of positive and negative IR z-scores
table(sign(df_IR_ref_z_plot$IR_z))


# ##############################
# # sample PCA 
# ################################

# plotPCA <- function(df_vals, list_vals, file_tag, meta_col, title ){

#   set.seed(123)

#   head(df_vals)
#   dim(df_vals)

#   df_vals_sigil <- df_vals %>%
#     filter(row.names(.) %in% list_vals)

#   head(df_vals)
#   dim(df_vals_sigil)

#   # remove high nan; in remaining rows replace with row median
#   df_vals_sigil_clean <- df_vals_sigil[which(rowMeans(!is.na(df_vals_sigil)) > 0.5), ]  %>%
#       mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = T), x))

#   print(dim(df_vals_sigil_clean))

#   df_vals_sigil_clean <- df_vals_sigil_clean[apply(df_vals_sigil_clean, 1, var) != 0, ]
#   print(dim(df_vals_sigil_clean))


#   prcomp.out <- prcomp(as.data.frame(t(df_vals_sigil_clean)),
#                     center = TRUE,
#                     scale. = TRUE)
    
#   # print(summary(prcomp.out))
#   # Variance explained by each PC
#   var_explained <- prcomp.out$sdev^2/sum(prcomp.out$sdev^2)

#   # Merge PCA results with metadata
#   df_PCA <- data.frame(x = prcomp.out$x[,1],  y = prcomp.out$x[,2])
#   rownames(df_PCA) <- colnames(df_vals_sigil_clean)
#   pca.out.merge = cbind(df_PCA, df_metadata)
#   print(dim(pca.out.merge ))


#   n <- length(unique(df_metadata[[meta_col]]))
#   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#   pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#   # Plot PCA all 
#   plt_pca <- ggplot(pca.out.merge, aes(x, y, color = get(meta_col), shape = data_source)) +
#       geom_point(size = 2.5) +
#       theme_classic() +
#       theme(legend.position="top",legend.title = element_blank(), legend.text = element_text(size = 10))+
#       # guides(color=guide_legend(nrow=3,byrow=TRUE)) +
#       scale_color_manual(values=pal) +
#       labs(title= paste0(title), sep = ' ', 
#           x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#           y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) 

#   # Save plot
#   ggsave(file.path(fig_output_path,
#                   paste0("PCA_",file_tag,".png")),
#           device = "png",
#           width = 5, height = 8,
#           dpi = 300)

#   # UMAP using PCA as input

#   # Run UMAP
#   umap.out <- umap(as.data.frame(prcomp.out$x), n_neighbors = 5, learning_rate = 0.5, init = "random", min_dist = 1, spread = 5)
#   umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
#   rownames(umap.out) <- colnames(df_vals_sigil_clean)

#   # Merge UMAP results with metadata
#   umap.out.merge = cbind(umap.out, df_metadata)

#   # Plot UMAP
#   plt_umap <- ggplot(umap.out.merge, aes(x, y, color = get(meta_col), shape = data_source)) +
#       geom_point(size = 2.5) +
#       theme_classic() +
#       theme(legend.position="top",legend.title = element_blank(), legend.text = element_text(size = 10))+
#       # guides(color=guide_legend(nrow=3,byrow=TRUE)) +
#       scale_color_manual(values=pal) +
#       labs(title= paste0(title), sep = ' ', 
#           x="UMAP1",
#           y="UMAP2") 

#   # Save plot
#   ggsave(file.path(fig_output_path,
#                   paste0("UMAP_",file_tag,".png")),
#           device = "png",
#           width = 5, height = 8,
#           dpi = 300)

#   return(list("pca"=plt_pca, "umap" = plt_umap))

# }

# ###################################################
# # GROUP LABEL SIGIL
# # Make PCAs and UMAPs by group label 
# gene_by_group <- plotPCA(df_exp,unique(df_gene_set$X), "Gene_sigil_by_group", "group_label", "Gene" )
# PS_by_group <- plotPCA(df_all_PS,unique(df_splice_set$event), "PS_sigil_by_group", "group_label", "Splice" )
# IR_by_group <- plotPCA(df_IR_table,unique(df_IR_set$event), "IR_sigil_by_group", "group_label", "Intron retention" )

# # Plot all by group label
# prow <- plot_grid(
#   gene_by_group$pca + theme(legend.position="none"),
#   PS_by_group$pca + theme(legend.position="none"),
#   IR_by_group$pca + theme(legend.position="none"),
#   align = 'vh',
#   # labels = c("A", "B", "C"),
#   hjust = -1,
#   nrow = 1
# )

# # Get legend from plot
# legend <- get_legend(
#   gene_by_group$pca + theme(legend.position="bottom")
# )

# # Add the legend to the row 
# sigil_by_group_pca <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .3))
# ggsave(file.path(fig_output_path,
#                 paste0("PCA_all_by_group.png")),
#         device = "png",
#         width = 15, height = 6,
#         dpi = 300, bg="white")

# #Plot UMAPS
# prow_umap <- plot_grid(
#   gene_by_group$umap + theme(legend.position="none"),
#   PS_by_group$umap + theme(legend.position="none"),
#   IR_by_group$umap + theme(legend.position="none"),
#   align = 'vh',
#   # labels = c("A", "B", "C"),
#   hjust = -1,
#   nrow = 1
# )
# sigil_by_group_umap <- plot_grid(prow_umap, legend, ncol = 1, rel_heights = c(1, .3))
# ggsave(file.path(fig_output_path,
#                 paste0("UMAP_all_by_group.png")),
#         device = "png",
#         width = 15, height = 6,
#         dpi = 300, bg="white")
# ###################################################
# # MAIN LABEL SIGIL
# # Make PCAs and UMAPs by main label 
# gene_by_main <- plotPCA(df_exp,unique(df_gene_set$X), "Gene_sigil_by_main", "main_label", "Gene" )
# PS_by_main <- plotPCA(df_all_PS,unique(df_splice_set$event), "PS_sigil_by_main", "main_label", "Splice" )
# IR_by_main <- plotPCA(df_IR_table,unique(df_IR_set$event), "IR_sigil_by_main", "main_label", "Intron retention" )


# # PCA Plot all by main label
# prow <- plot_grid(
#   gene_by_main$pca + theme(legend.position="none"),
#   PS_by_main$pca + theme(legend.position="none"),
#   IR_by_main$pca + theme(legend.position="none"),
#   align = 'vh',
#   # labels = c("A", "B", "C"),
#   hjust = -1,
#   nrow = 1
# )

# # Get legend from plot
# legend <- get_legend(
#   gene_by_main$pca + theme(legend.position="bottom")
# )

# # PCA Add the legend to the row 
# sigil_by_main_pca <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .3))
# ggsave(file.path(fig_output_path,
#                 paste0("PCA_all_by_main.png")),
#         device = "png",
#         width = 15, height = 6,
#         dpi = 300, bg="white")
# #Plot UMAPS
# prow_umap <- plot_grid(
#   gene_by_main$umap + theme(legend.position="none"),
#   PS_by_main$umap + theme(legend.position="none"),
#   IR_by_main$umap + theme(legend.position="none"),
#   align = 'vh',
#   # labels = c("A", "B", "C"),
#   hjust = -1,
#   nrow = 1
# )
# sigil_by_main_umap <- plot_grid(prow_umap, legend, ncol = 1, rel_heights = c(1, .3))
# ggsave(file.path(fig_output_path,
#                 paste0("UMAP_all_by_main.png")),
#         device = "png",
#         width = 15, height = 6,
#         dpi = 300, bg="white")
# ###################################################
# # GROUP LABEL ALL JUNCTIONS/GENES
# # Make PCAs by group label 
# gene_by_group_all <- plotPCA(df_exp,rownames(df_exp), "Gene_sigil_by_group", "group_label", "Gene" )
# PS_by_group_all <- plotPCA(df_all_PS,rownames(df_all_PS), "PS_sigil_by_group", "group_label", "Splice" )
# IR_by_group_all <- plotPCA(df_IR_table,rownames(df_IR_table), "IR_sigil_by_group", "group_label", "Intron retention" )

# # Plot all by group label
# prow <- plot_grid(
#   gene_by_group_all$pca + theme(legend.position="none"),
#   PS_by_group_all$pca + theme(legend.position="none"),
#   IR_by_group_all$pca + theme(legend.position="none"),
#   align = 'vh',
#   # labels = c("A", "B", "C"),
#   hjust = -1,
#   nrow = 1
# )

# # Get legend from plot
# legend <- get_legend(
#   gene_by_group_all$pca + theme(legend.position="bottom")
# )

# # Add the legend to the row 
# all_genes_juncs_by_group_pca <-plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .3))
# ggsave(file.path(fig_output_path,
#                 paste0("PCA_all_by_group_all_juncs_or_genes.png")),
#         device = "png",
#         width = 15, height = 6,
#         dpi = 300, bg="white")

# prow_umap <- plot_grid(
#   gene_by_group_all$umap + theme(legend.position="none"),
#   PS_by_group_all$umap + theme(legend.position="none"),
#   IR_by_group_all$umap + theme(legend.position="none"),
#   align = 'vh',
#   # labels = c("A", "B", "C"),
#   hjust = -1,
#   nrow = 1
# )
# all_genes_juncs_by_group_umap <- plot_grid(prow_umap, legend, ncol = 1, rel_heights = c(1, .3))
# ggsave(file.path(fig_output_path,
#                 paste0("UMAP_all_by_group_all_juncs_or_genes.png")),
#         device = "png",
#         width = 15, height = 6,
#         dpi = 300, bg="white")
# ###################################################
# # Combining figure rows

# # Combine all 3 PCA fig rows
# plot_grid(all_genes_juncs_by_group_pca,sigil_by_group_pca, sigil_by_main_pca, ncol = 1)
#         # labels = c("All genes or junctions", "Sigil genes or junctions", "Sigil genes or junctions"), label_x = 0)
#         # labels = c("All genes or junctions", "Sigil genes or junctions", "Sigil genes or junctions")

# ggsave(file.path(fig_output_path,
#                 paste0("PCA_all_figrows.png")),
#         device = "png",
#         width = 15, height = 22,
#         dpi = 300, bg="white")


# # Combine all 3 UMAP fig rows
# plot_grid(all_genes_juncs_by_group_umap,sigil_by_group_umap, sigil_by_main_umap, ncol = 1)
# ggsave(file.path(fig_output_path,
#                 paste0("UMAP_all_figrows.png")),
#         device = "png",
#         width = 15, height = 22,
#         dpi = 300, bg="white")