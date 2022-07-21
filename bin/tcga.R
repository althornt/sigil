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
library(purrr)
library(ComplexHeatmap)
library(UpSetR)
library(cowplot)
library(circlize)

# library(enrichR)
# library(GenomicRanges)
# library(valr)
# library(rGREAT)

sigil_out_path = "/mnt/results_sigil_combine/sigil_results_SongChoi_newlabels_20220614/"
fig_output_path = "/mnt/tcga/outputs/"
if (!dir.exists(fig_output_path)){
    dir.create(fig_output_path,
    recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"), stringsAsFactors=FALSE)
# print(df_metadata)

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

# df_IR_set <- read.csv(file = paste0(sigil_out_path,
#                           "combine_mesa_out/splice_set_IR/df_splice_sets.csv"))
# print(head(df_IR_set))

# df_gene_set <- read.csv(file = paste0(sigil_out_path,
#                           "combine_gene_out/gene_sets/df_gene_sets.csv"))
# print(head(df_gene_set))

# # # Read in the huge file and save version with only splice set junctions 
# df_LUAD <- read.table(file = "/mnt/tcga/2022.07.06.luad_allPS.tsv",
#                           sep="\t", header = TRUE, row.names=1, check.names = FALSE) 
                        
# #reformat events/rownames to drop "chr"
# rownames(df_LUAD) <- str_remove(rownames(df_LUAD), "chr")

# # print(head(df_LUAD))
# print(dim(df_LUAD))
# # print(head(unique(df_splice_set$event)))

# df_LUAD_spliceset <- df_LUAD %>%
#     filter(row.names(.) %in% unique(df_splice_set$event))

# # print(dim(df_LUAD_spliceset))

# write.table(df_LUAD_spliceset,
#           na="nan", row.names = TRUE, quote=FALSE,
#            sep = "\t",
#           file ="/mnt/tcga/2022.07.06.luad_allPS_spliceset.tsv")

df_LUAD_spliceset <- read.table(file = "/mnt/tcga/2022.07.06.luad_allPS_spliceset.tsv",
                          sep="\t", header = TRUE, row.names=1, check.names = FALSE) %>%
                          filter( rowSums(is.na(.)) != ncol(.)) #drop rows where all are NA

print(head(colnames(df_LUAD_spliceset)))


# print(head(df_LUAD_spliceset))
cat("Dim of LUAD PS only splice set ")
print(dim(df_LUAD_spliceset))


# Read in and merge tcga metadata

# File to map the sample names which are the bam ids to clincial 
df_bam_id <- read.table(file = "/mnt/tcga/gdc_sample_sheet.2022-07-08.tsv",
                          sep="\t", header = TRUE)

print(head(df_bam_id))

df_bam_id$File.Name <- str_remove(df_bam_id$File.Name, ".rna_seq.genomic.gdc_realn.bam")
print(head(df_bam_id))

# make mesa manifests
# tumor_samples <- df_bam_id %>%
#     filter(Sample.Type %in% c("Primary Tumor", "Recurrent Tumor")) %>%
#     select(File.Name)
# print(length(tumor_samples))
# # quit()
# write.table(x = tumor_samples,row.names = FALSE, quote=FALSE,col.names=FALSE,
#         file = paste0(paste0(fig_output_path,"manifest_LUAD_tumor.txt")))

# normal_samples <- df_bam_id %>%
#     filter(Sample.Type == "Solid Tissue Normal") %>%
#     select(File.Name)
# print(length(normal_samples))

# write.table(x = normal_samples,row.names = FALSE, quote=FALSE,col.names=FALSE,
#         file = paste0(paste0(fig_output_path,"manifest_LUAD_normal.txt")))


# Read in mesa normal vs tumor results 
df_mesa_tumor_vs_norm <- read.table(
        file = "/mnt/tcga/outputs/LUAD_normal_vs_tumor_mesa_comp.txt",
        sep="\t", header = TRUE) %>%
    arrange(corrected)

df_mesa_tumor_vs_norm$event <- str_remove(df_mesa_tumor_vs_norm$event, "chr")

# Rename the tcga associated columns
df_mesa_tumor_vs_norm_sigil <- df_mesa_tumor_vs_norm %>%
    filter(event %in% unique(df_splice_set$event))  %>%
    filter((corrected < .05) & (abs(delta) > .1)) %>%
    rename(normal_mean = mean1) %>%
    rename(cancer_mean = mean2) %>%
    rename(norm_vs_cancer_delta = delta ) %>%
    rename(norm_vs_cancer_corrected_pval = corrected) %>%
    column_to_rownames("event") %>%
    select(normal_mean, cancer_mean,norm_vs_cancer_delta, norm_vs_cancer_corrected_pval ) 

print(dim(df_mesa_tumor_vs_norm_sigil))

# Merge tcga info into splice set
df_splice_set <- merge(df_splice_set, df_mesa_tumor_vs_norm_sigil, by.x = "event", by.y = 0,sort=FALSE, all= TRUE) 
df_splice_set <- df_splice_set %>%
    mutate(delta_match= ifelse(sign(delta) == sign(norm_vs_cancer_delta), "match", "no_match" )) 

# Make df with row per event listing sets 
df_splice_set_byevent <- df_splice_set %>%
    group_by(event) %>%
    summarise(ls_sets = toString(set)) %>%
    as.data.frame() %>%
    column_to_rownames("event") 

print("splcie set by event")
dim(df_splice_set_byevent)

# Merge tcga and splice set 
df_splice_set_byevent <- merge(df_splice_set_byevent, df_mesa_tumor_vs_norm_sigil, by=0, all=FALSE)
print(dim(df_splice_set_byevent))

df_splice_set_byevent <- merge(df_splice_set_byevent, df_splice_set, by=0, all=FALSE)

print("merge")
dim(df_splice_set_byevent)
head(df_splice_set_byevent)

cat("Number of events: \n")
print(length(unique(df_splice_set_byevent$event)))
cat("Number of genes: \n")
print(length(unique(df_splice_set_byevent$overlapping)))

# for (i in sort(unique(df_splice_set_byevent$overlapping))){
#     cat("\n")
#     cat( i )
# }

###########################
# Merging metadata 
###############################
# df_LUAD_dcc <- read.table(file = "/mnt/tcga/dcc_sample_data_scrape_all.tsv",
#                           sep="\t", header = TRUE) %>%
#                     filter(cases.0.project.project_id == "TCGA-LUAD") 
df_tcga_clinical <- read.table(file = "/mnt/tcga/luad_tcga_pan_can_atlas_2018_clinical_data.tsv",
                          sep="\t", header = TRUE) %>%
                          filter(Patient.ID %in% df_bam_id$Case.ID)

df_meta <- merge(df_bam_id, df_tcga_clinical, by.x = "Case.ID", by.y = "Patient.ID",sort=FALSE, all= TRUE) 

print(length(df_bam_id$Case.ID))
print(length(unique(df_bam_id$Case.ID)))

print(length(unique(df_meta$Case.ID)))
print(length(df_meta$Case.ID))
print(head(df_meta$Case.ID))
# quit()
cat("Dim of bam meta")
print(dim(df_bam_id))

cat("Dim of clinical meta")
print(dim(df_tcga_clinical))

cat("Dim merged metadata ")
dim(df_meta)

cat("Intersection of samples in metadata and PS table")
print(length(intersect(colnames(df_LUAD_spliceset), df_meta$File.Name)))

cat("Intersection of samples in bam metadata and PS table")
print(length(intersect(colnames(df_LUAD_spliceset),df_bam_id$File.Name)))

# df_meta %>%
#     filter(Sample.Type.x == "Solid Tissue Normal") %>%
#     select(c("Sample.Type.x", "Case.ID")) %>%

colnames(df_meta)

df_meta <- df_meta %>%
    filter(File.Name %in% colnames(df_LUAD_spliceset)) %>%
    # distinct() %>%
    # select(c("Tumor.Type", "File.Name", "MSI.MANTIS.Score", "Mutation.Count", 
    #             "Person.Neoplasm.Cancer.Status", "Ragnum.Hypoxia.Score", 
    #             "Sex", "Sample.Type.x", "Cancer.Type", "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code")) %>%
    column_to_rownames("File.Name") %>%
    rename(LUAD.Sample.Type = Sample.Type.x) %>%
    arrange(rownames(.))

sort(colnames(df_meta))
head(df_meta)
# quit()
dim(df_meta)
df_bam_id <- df_bam_id %>%
    filter(File.Name %in% colnames(df_LUAD_spliceset)) %>%
    column_to_rownames("File.Name") %>%
    arrange(rownames(.))
head(df_bam_id)

# GSVA
df_enr_all_sets <- read.csv(file = "/mnt/tcga/GSVA_TCGA_LUAD_output.tsv", stringsAsFactors = TRUE, check.names = FALSE) %>%
    tibble::column_to_rownames("name") %>%
        select( intersect(colnames(.),rownames(df_meta))) %>%
        select(sort(names(.)))

head(colnames(df_enr_all_sets))

df_enr_group_sets <- read.csv(file = "/mnt/tcga/GSVA_TCGA_LUAD_output_group_set.tsv", stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames("name") %>%
    select( intersect(colnames(.),rownames(df_meta))) %>%
    select(sort(names(.)))

head(colnames(df_enr_all_sets))
head(rownames(df_meta))
################################################
# Counts matching sign of delta 
################################################
total_potential <- df_splice_set %>% 
    filter(!is.na(norm_vs_cancer_delta)) %>%
    nrow()

cat("\n number of potential LUAD junction set matches: \n")
cat(total_potential, "\n")

matches <- df_splice_set    %>%
    filter(delta_match == "match") %>%
    select(event, delta_match, set) %>%
    arrange(set)
    # select( norm_vs_cancer_delta, delta)

cat("\n number of LUAD junction set matches with same delta: \n")
cat(nrow(matches), "\n")


cat("\n percent of LUAD junction set matches with same delta: \n")
cat(100*(nrow(matches)/total_potential), "\n")


cat("\n count unique junctions that match a set delta: \n")
cat(length(unique(matches$event)), "\n")


################################################
# volcano plot 
################################################

# head(df_mesa_tumor_vs_norm)

# print(length(unique(df_splice_set$event)))
# print(length(unique(matches$event)))

# df_mesa_tumor_vs_norm <- df_mesa_tumor_vs_norm %>%
#     mutate(sigil = ifelse(event %in% unlist(df_splice_set$event), 
#         "Junction in a sigil splice set", "Junction not in a sigil splice set")) %>%
#     mutate(sigil_and_delta_match = ifelse(event %in% matches$event, 
#         "Junction in a sigil splice set with matching delta direction", "Junction not in a sigil splice set with matching delta direction"))

# print(length(df_mesa_tumor_vs_norm$sigil))
# print(length(df_mesa_tumor_vs_norm$sigil_and_delta_match))

# # Plot coloring all sigil 
# p <- ggplot(data=df_mesa_tumor_vs_norm, aes(x=delta, y=-log10(corrected), color=sigil )) +
#         geom_point(alpha = .2  ) +
#         # theme_minimal() + geom_vline(xintercept=c(-0.2, 0.2), col="red") +
#         # geom_hline(yintercept=-log10(0.05), col="red") +
#         xlim(-1.0, 1.0) +
#         theme_classic()+
#         # +
#         # geom_text(
#         # label= df_css$gglabel,
#         # vjust="inward",hjust="inward",
#         # # nudge_x = 0.05, nudge_y = 0.05,
#         # check_overlap =F, col = "darkgreen", size = 2
#         # ) +
#         scale_color_manual(name = "",
#         values = c("Junction in a sigil splice set" = "red",
#                     "Junction not in a sigil splice set" = "black"),
#         labels = c("Junction in a sigil splice set",
#                     "Junction not in a sigil splice set"))+
#         theme(legend.position="bottom") +
# #   theme(
# #     legend.position="bottom", 
# #     legend.justification = "left", 
# #     legend.margin = margin(0, 0, 0, 0),
# #     legend.spacing.x = unit(0, "pt")
# #   )+
#           guides(color = guide_legend(nrow = 3)) +
#         ylab("-log10(corrected p-value")


# ggsave(file.path(fig_output_path,
#                 paste("LUAD_volcano_all_sigil_juncs.png", sep = '.')),
#         device = "png",
#         width = 5, height = 4,
#         dpi = 300)

# # Plot coloring all sigil sig and matched
# p <- ggplot(data=df_mesa_tumor_vs_norm, aes(x=delta, y=-log10(corrected), color=sigil_and_delta_match )) +
#         geom_point(alpha = .2  ) +
#         # theme_minimal() + geom_vline(xintercept=c(-0.2, 0.2), col="red") +
#         # geom_hline(yintercept=-log10(0.05), col="red") +
#         xlim(-1.0, 1.0) +
#         theme_classic()+
#         # +
#         # geom_text(
#         # label= df_css$gglabel,
#         # vjust="inward",hjust="inward",
#         # # nudge_x = 0.05, nudge_y = 0.05,
#         # check_overlap =F, col = "darkgreen", size = 2
#         # ) +
#         scale_color_manual(name = "",
#         values = c("Junction in a sigil splice set with matching delta direction" = "red",
#                     "Junction not in a sigil splice set with matching delta direction" = "black"),
#         labels = c("Junction in a sigil splice set with matching delta direction",
#                     "Junction not in a sigil splice set with matching delta direction"))+
#         theme(legend.position="bottom") +
#         guides(color = guide_legend(nrow = 3)) +
#         ylab("-log10(corrected p-value)") 

# ggsave(file.path(fig_output_path,
#                 paste("LUAD_volcano_sigil_sig_delta_match.png", sep = '.')),
#         device = "png",
#         width = 5, height = 4,
#         dpi = 300)


# ################################################
# # Heatmap sigil set
# ################################################
# df_LUAD_spliceset_high_var <- df_LUAD_spliceset %>%
#     select( intersect(colnames(.),rownames(df_meta))) %>%
#     # filter(rownames(.) %in% rownames(df_mesa_tumor_vs_norm_sigil))%>%
#     rowwise() %>%
#     mutate(variance = c_across(everything()) %>% var()) 

# rownames(df_LUAD_spliceset_high_var) <- rownames(df_LUAD_spliceset)
# df_LUAD_spliceset_high_var <- df_LUAD_spliceset_high_var %>%
#      ungroup() %>%
#      arrange(desc(variance)) %>%
#      head(300) %>%
#     as.data.frame() %>%
#     select(!variance) %>%
#     select(sort(names(.)))

# # print(head(df_LUAD_spliceset_high_var))
# print(dim(df_LUAD_spliceset_high_var))
# # quit()
# df_LUAD_spliceset_scaled = t(scale(t(df_LUAD_spliceset_high_var)))

# # Check order match
# head(colnames(df_LUAD_spliceset_scaled))
# head(rownames(df_meta))


# ha <- HeatmapAnnotation(
#     df = df_meta %>% select(c("LUAD.Sample.Type")) , 
#     col = list(LUAD.Sample.Type = c("Primary Tumor" ="#F8766D", 
#                             "Solid Tissue Normal"= "#619CFF",
#                              "Recurrent Tumor"="#00BA38"    )
#                 # Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code = c(
#                 #                     "NA" = "gray" ,
#                 #                     "STAGE IB" = "black" ,
#                 #                     "STAGE I" = "black" ,
#                 #                     "STAGE IA" = "black" ,
#                 #                     "STAGE II" = "yellow" ,
#                 #                     "STAGE IIA" = "yellow" ,
#                 #                     "STAGE IIB" = "yellow" ,
#                 #                     "STAGE IIIA" = "orange" ,
#                 #                     "STAGE IIIB" = "orange" ,
#                 #                     "STAGE IV" = "red" 
            
#                 )            )

# ht <- ComplexHeatmap::Heatmap(df_LUAD_spliceset_scaled,
#                     name = "Z-Score PS",
#                 # cluster_columns = FALSE,
#                                 # row_title = "UP", row_title_rot = 0,
#                                 # column_order =order(colnames(df_enr_median_heat_UP)),
#                                 # row_order = order(rownames(df_enr_median_heat_UP)), 
#                                 show_row_names= FALSE,
#                                 show_column_names = FALSE,
#                                 # row_names_gp = grid::gpar(fontsize =8),
#                                 # column_names_gp = grid::gpar(fontsize =7),
#                                 # show_heatmap_legend = FALSE
#                                 top_annotation=ha,
#                                 # show_row_dend = TRUE,
#                                 # heatmap_legend_param = list(
#                                 # legend_direction = "horizontal", 
#                                 # legend_height = unit(1, "cm"),
#                                 # legend_gp = gpar(fontsize = 5)
#                                 )
# # combined gene and splice heatmaps
# png(file=paste0(fig_output_path,"LUAD_spliceset_heatmap.png"),
#     width = 50,
#     height    = 20,
#     units     = "cm",
#     res       = 1200)

# draw(ht, legend_grouping = "original",heatmap_legend_side="bottom", annotation_legend_side="bottom")
# dev.off()

###################################################
# PCA tumor vs normal + sigil set
####################################################

df_LUAD_spliceset_clean <- df_LUAD_spliceset[which(rowMeans(!is.na(df_LUAD_spliceset)) > 0.5), ]  %>%
    select(intersect(colnames(df_LUAD_spliceset),rownames(df_meta))) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = T), x)) %>%
    select(sort(names(.))) 

print(dim(df_LUAD_spliceset_clean))

prcomp.out <- prcomp(as.data.frame(t(df_LUAD_spliceset_clean)),
                   center = TRUE,
                   scale. = TRUE)

var_explained <- prcomp.out$sdev^2/sum(prcomp.out$sdev^2)

# print(head(prcomp.out ))
print(dim(prcomp.out ))

# Merge PCA results with metadata
df_PCA <- data.frame(x = prcomp.out$x[,1],  y = prcomp.out$x[,2])
rownames(df_PCA) <- colnames(df_LUAD_spliceset_clean)
pca.out.merge = cbind(df_PCA, df_meta)
print(dim(pca.out.merge ))

# # Make color palette
# n <- length(unique(df_bam_id[["Sample.Type"]]))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Plot PCA all s
plt <- ggplot(pca.out.merge, aes(x, y, color = LUAD.Sample.Type)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="top",legend.title = element_blank()) +
    # scale_color_manual(values=pal) +
    labs(title= "", sep = ' ', x="PC1", y="PC2"
        #   x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
        #   y=paste0("PC2: ",round(var_explained[2]*100,1),"%")
          ) 

# Save plot
ggsave(file.path(fig_output_path,
                paste("LUAD_spliceset_PCA.png", sep = '.')),
        device = "png",
        width = 5, height = 4,
        dpi = 300)

# ############################################################################
# ################################################
# # Heatmap tumor vs normal + sigil set + matched delta
# ################################################
# df_LUAD_sig_diff_and_sigil <- df_LUAD_spliceset %>%
#     select( intersect(colnames(.),rownames(df_meta))) %>%
#     # filter(rownames(.) %in% rownames(df_mesa_tumor_vs_norm_sigil))%>%
#     filter(rownames(.) %in% matches$event)%>%
#     select(sort(names(.))) 

# df_LUAD_sig_diff_and_sigil_scaled = t(scale(t(df_LUAD_sig_diff_and_sigil)))

# ha <- HeatmapAnnotation(
#     df = df_meta %>% select(c("LUAD.Sample.Type")), 
#     col = list(LUAD.Sample.Type = c("Primary Tumor" ="#F8766D", 
#                             "Solid Tissue Normal"= "#619CFF",
#                              "Recurrent Tumor"="#00BA38"    )
#                 # Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code = c(
#                 #                     "NA" = "gray" ,
#                 #                     "STAGE IB" = "black" ,
#                 #                     "STAGE I" = "black" ,
#                 #                     "STAGE IA" = "black" ,
#                 #                     "STAGE II" = "yellow" ,
#                 #                     "STAGE IIA" = "yellow" ,
#                 #                     "STAGE IIB" = "yellow" ,
#                 #                     "STAGE IIIA" = "orange" ,
#                 #                     "STAGE IIIB" = "orange" ,
#                 #                     "STAGE IV" = "red" 
            
#                 )            )

# ht <- ComplexHeatmap::Heatmap(df_LUAD_sig_diff_and_sigil_scaled,
#                     name = "Z-Score PS",
#                 # cluster_columns = FALSE,
#                                 # row_title = "UP", row_title_rot = 0,
#                                 # column_order =order(colnames(df_enr_median_heat_UP)),
#                                 # row_order = order(rownames(df_enr_median_heat_UP)), 
#                                 show_row_names= FALSE,
#                                 show_column_names = FALSE,
#                                 # row_names_gp = grid::gpar(fontsize =8),
#                                 # column_names_gp = grid::gpar(fontsize =7),
#                                 # show_heatmap_legend = FALSE
#                                 top_annotation=ha,
#                                 # show_row_dend = TRUE,
#                                 heatmap_legend_param = list(
#                                 # legend_direction = "horizontal", 
#                                 legend_height = unit(10, "cm"),
#                                 legend_gp = gpar(fontsize = 20)
#                                 ))
# # combined gene and splice heatmaps
# png(file=paste0(fig_output_path,"LUAD_sig_diff_sigil_heatmap_scaled.png"),
#     width = 60,
#     height    = 20,
#     units     = "cm",
#     res       = 1200)

# draw(ht, legend_grouping = "original",heatmap_legend_side="right", annotation_legend_side="right")

# ht <- ComplexHeatmap::Heatmap(df_LUAD_sig_diff_and_sigil,
#                     name = "PS",
#                       col = circlize::colorRamp2(c(0, 1.0), c("white", "red")),

#                 # cluster_columns = FALSE,
#                                 # row_title = "UP", row_title_rot = 0,
#                                 # column_order =order(colnames(df_enr_median_heat_UP)),
#                                 # row_order = order(rownames(df_enr_median_heat_UP)), 
#                                 show_row_names= FALSE,
#                                 show_column_names = FALSE,
#                                 # row_names_gp = grid::gpar(fontsize =8),
#                                 # column_names_gp = grid::gpar(fontsize =7),
#                                 # show_heatmap_legend = FALSE
#                                 top_annotation=ha,
#                                 # show_row_dend = TRUE,
#                                 heatmap_legend_param = list(
#                                 # legend_direction = "horizontal", 
#                                 legend_height = unit(2, "cm"),
#                                 # legend_gp = gpar(fontsize = 5), 
#                                 legend_gp = gpar(fontsize = 14)
#                                 ))
# # combined gene and splice heatmaps
# png(file=paste0(fig_output_path,"LUAD_sig_diff_sigil_heatmap_unscaled.png"),
#     width = 60,
#     height    = 20,
#     units     = "cm",
#     res       = 1200)

# draw(ht, legend_grouping = "original",heatmap_legend_side="right", annotation_legend_side="right")
# dev.off()

# ###################################################
# # PCA tumor vs normal + sigil set + matched delta
# ####################################################

# df_LUAD_sig_diff_and_sigil_clean <- df_LUAD_sig_diff_and_sigil[which(rowMeans(!is.na(df_LUAD_sig_diff_and_sigil)) > 0.5), ]  %>%
#     select(intersect(colnames(df_LUAD_sig_diff_and_sigil),rownames(df_meta))) %>%
#     mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = T), x))

# print(dim(df_LUAD_sig_diff_and_sigil))

# prcomp.out <- prcomp(as.data.frame(t(df_LUAD_sig_diff_and_sigil_clean)),
#                    center = TRUE,
#                    scale. = TRUE)$x
  
# # print(head(prcomp.out ))
# print(dim(prcomp.out ))

# # Merge PCA results with metadata
# df_PCA <- data.frame(x = prcomp.out[,1],  y = prcomp.out[,2])
# rownames(df_PCA) <- colnames(df_LUAD_sig_diff_and_sigil_clean)
# pca.out.merge = cbind(df_PCA, df_meta)
# print(dim(pca.out.merge ))

# # # Make color palette
# # n <- length(unique(df_bam_id[["Sample.Type"]]))
# # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# # pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# # Plot PCA all s
# plt <- ggplot(pca.out.merge, aes(x, y, color = LUAD.Sample.Type)) +
#     geom_point(size = 2) +
#     theme_classic() +
#     theme(legend.position="top",legend.title = element_blank()) +
#     # scale_color_manual(values=pal) +
#     labs(title= "", sep = ' ')

# # Save plot
# ggsave(file.path(fig_output_path,
#                 paste("LUAD_sig_diff_sigil_PCA.png", sep = '.')),
#         device = "png",
#         width = 5, height = 4,
#         dpi = 300)

#################################
# Mapping sets and cell types
#################################
Bcell_sets <- list(
    "B_cells_group_DN",
    "B_cells_group_UP",
    "B_cells_memory_main_DN",
    "B_cells_memory_main_UP",
    "B_cells_memory_within_DN",
    "B_cells_memory_within_UP",
    "B_cells_naive_main_DN",
    "B_cells_naive_main_UP",
    "B_cells_naive_within_DN",
    "B_cells_naive_within_UP"
)

Tcell_sets <- list(
   "T_Cells_CD4:+_main_DN",
   "T_Cells_CD4:+_main_UP",
   "T_Cells_CD4:+_within_DN",
   "T_Cells_CD4:+_within_UP",
   "T_Cells_CD8:+_main_DN",
   "T_Cells_CD8:+_main_UP",
   "T_Cells_CD8:+_within_DN",
   "T_Cells_CD8:+_within_UP",
   "T_cells_group_DN",
   "T_cells_group_UP"
)
Macrophage_sets <- list(
    "Macrophage_LPS-18h_main_DN",
    "Macrophage_LPS-18h_main_UP",
    "Macrophage_LPS-18h_within_DN",
    "Macrophage_LPS-18h_within_UP",
    "Macrophage_NT_main_DN",
    "Macrophage_NT_main_UP",
    "Macrophage_NT_within_DN",
    "Macrophage_NT_within_UP",
    "Macrophage_Pam3CSK4-18h_main_DN",
    "Macrophage_Pam3CSK4-18h_main_UP",
    "Macrophage_Pam3CSK4-18h_within_DN",
    "Macrophage_Pam3CSK4-18h_within_UP",
    "Macrophage_R837-18h_main_DN",
    "Macrophage_R837-18h_main_UP",
    "Macrophage_R837-18h_within_DN",
    "Macrophage_R837-18h_within_UP",
    "Macrophage_R848-18h_main_DN",
    "Macrophage_R848-18h_main_UP",
    "Macrophage_R848-18h_within_DN",
    "Macrophage_R848-18h_within_UP",
    "Macrophages_group_DN",
    "Macrophages_group_UP"
)
DC_sets <- list(
 "DC_Myeloid_CD123+_main_DN",
 "DC_Myeloid_CD123+_main_UP",
 "DC_Myeloid_CD123+_within_DN",
 "DC_Myeloid_CD123+_within_UP",
 "DC_Myeloid_main_DN",
 "DC_Myeloid_main_UP",
 "DC_Myeloid_within_DN",
 "DC_Myeloid_within_UP",
 "DC_Plasmacytoid_main_DN",
 "DC_Plasmacytoid_main_UP",
 "DC_Plasmacytoid_within_DN",
 "DC_Plasmacytoid_within_UP",
 "Dendritic_LPS-18h_main_DN",
 "Dendritic_LPS-18h_main_UP",
 "Dendritic_LPS-18h_within_DN",
 "Dendritic_LPS-18h_within_UP",
 "Dendritic_NT_main_DN",
 "Dendritic_NT_main_UP",
 "Dendritic_NT_within_DN",
 "Dendritic_NT_within_UP",
 "Dendritic_R837-18h_main_DN",
 "Dendritic_R837-18h_main_UP",
 "Dendritic_R837-18h_within_DN",
 "Dendritic_R837-18h_within_UP",
 "Dendritic_R848-18h_main_DN",
 "Dendritic_R848-18h_main_UP",
 "Dendritic_R848-18h_within_DN",
 "Dendritic_R848-18h_within_UP",
 "Dendritic_cells_group_DN",
 "Dendritic_cells_group_UP"
)

Eos_sets <- list(
    "Eosinophils_group_DN",
    "Eosinophils_group_UP"
)


Monocyte_sets <- list(
    "Monocyte_Choi_main_DN",
    "Monocyte_Choi_main_UP",
    "Monocyte_Choi_within_DN",
    "Monocyte_Choi_within_UP",
    "Monocyte_LPS-18h_main_DN",
    "Monocyte_LPS-18h_main_UP",
    "Monocyte_LPS-18h_within_DN",
    "Monocyte_LPS-18h_within_UP",
    "Monocyte_NT_main_DN",
    "Monocyte_NT_main_UP",
    "Monocyte_NT_within_DN",
    "Monocyte_NT_within_UP",
    "Monocyte_Pam3CSK4-18h_main_DN",
    "Monocyte_Pam3CSK4-18h_main_UP",
    "Monocyte_Pam3CSK4-18h_within_DN",
    "Monocyte_Pam3CSK4-18h_within_UP",
    "Monocyte_R837-18h_main_DN",
    "Monocyte_R837-18h_main_UP",
    "Monocyte_R837-18h_within_DN",
    "Monocyte_R837-18h_within_UP",
    "Monocyte_R848-18h_main_DN",
    "Monocyte_R848-18h_main_UP",
    "Monocyte_R848-18h_within_DN",
    "Monocyte_R848-18h_within_UP",
    "Monocytes_main_DN",
    "Monocytes_main_UP",
    "Monocytes_within_DN",
    "Monocytes_within_UP",
    "Monocytes_non-classical_main_DN",
    "Monocytes_non-classical_main_UP",
    "Monocytes_non-classical_within_DN",
    "Monocytes_non-classical_within_UP",
    "Monocytes_group_DN",
    "Monocytes_group_UP"
)

NK_sets <- list(
    "NK_cells_group_DN",
    "NK_cells_group_UP"
)

Neutrophil_sets <- list(
    "Neutrophils_group_DN",
    "Neutrophils_group_UP"
)

ls_map_set_type <- list(
  "Tcells" = list(Tcell_sets),
  "Bcells" = list(Bcell_sets),
  "Macrophages" = list(Macrophage_sets),
  "Monocytes" = list( Monocyte_sets),
  "NKcells" = list( NK_sets),
  "DC" = list(DC_sets),
  "Neutrophils"    =  list( Neutrophil_sets),
    "Eosinophils"    =  list( Eos_sets)

)
################################################
# Counts of sets in the sigil LUAD junctions 
################################################

# print(length(unique(unlist(ls_map_set_type))))


# matches_count <- matches %>%
#     count(set) %>%
#     arrange(desc(n)) %>%
#     head(20)

# ggplot(matches_count, aes(x=reorder(set,desc(n)), y=n)) +
#     geom_bar(stat = 'identity', position = "dodge") +
#     theme_classic() +
#     labs(x ="Set", y = "Number of events") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         text = element_text(size = 15)) 

#             # scale_fill_manual("name", values=c("#999999", "black"))+
#     # theme(legend.title= element_blank(), axis.title.x=element_blank()) +
#     # labs(title = paste0("Fishers test p-value = ", signif(res_fishers$p.value, 4)))

# ggsave(paste0(fig_output_path,"barplot_tcga_delta_matched_set_counts.png"), width=20, height=20, unit="cm", dpi = 300)

# quit()
# print(length(unique(df_splice_set$set)))
# matches$set_cat <- "NA"

# matches_count_by_cat <- matches %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[1]]), "T cells", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[2]]), "B cells", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[3]]), "Macrophages", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[4]]), "Monocytes", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[5]]), "NK cells", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[6]]), "Dendritic Cells", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[7]]), "Neutrophils", set_cat)) %>%
#     mutate(set_cat = ifelse(set %in% unlist(ls_map_set_type[[8]]), "Eosinophils", set_cat)) %>%
#     count(set_cat)

# head(matches_count_by_cat)

# ggplot(matches_count_by_cat, aes(x=set_cat, y=n)) +
#     geom_bar(stat = 'identity', position = "dodge") +
#     theme_classic() +
#     labs(x ="Set category", y = "Number of events") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         text = element_text(size = 15)) 

#             # scale_fill_manual("name", values=c("#999999", "black"))+
#     # theme(legend.title= element_blank(), axis.title.x=element_blank()) +
#     # labs(title = paste0("Fishers test p-value = ", signif(res_fishers$p.value, 4)))

# ggsave(paste0(fig_output_path,"barplot_tcga_delta_matched_set_counts_by_cat.png"), width=20, height=10, unit="cm", dpi = 300)



############################
# GSVA heatmap
###########################
# # check sample orders match
# head(colnames(df_enr_all_sets))
# head(rownames(df_meta))

# df_enr_all_sets_scaled = t(scale(t(df_enr_all_sets)))
# df_enr_group_sets_scaled = t(scale(t(df_enr_group_sets)))

# dim(df_enr_all_sets_scaled)
# dim(df_meta)

# ha <- HeatmapAnnotation(
#     # df = df_meta %>% select(c("LUAD.Sample.Type", "Tumor.Type")) , 
#     df = df_meta %>% select(c("LUAD.Sample.Type")) , 
#     col = list(LUAD.Sample.Type = c("Primary Tumor" ="#F8766D", 
#                             "Solid Tissue Normal"= "#619CFF",
#                              "Recurrent Tumor"="#00BA38"    )
#                 # Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code = c(
#                 #                     "NA" = "gray" ,
#                 #                     "STAGE IB" = "black" ,
#                 #                     "STAGE I" = "black" ,
#                 #                     "STAGE IA" = "black" ,
#                 #                     "STAGE II" = "yellow" ,
#                 #                     "STAGE IIA" = "yellow" ,
#                 #                     "STAGE IIB" = "yellow" ,
#                 #                     "STAGE IIIA" = "orange" ,
#                 #                     "STAGE IIIB" = "orange" ,
#                 #                     "STAGE IV" = "red" 
            
#                 )            )

# ht <- ComplexHeatmap::Heatmap(df_enr_all_sets_scaled,
#                     name = "Z-Score enrichment score",
#                 # cluster_columns = FALSE,
#                                 # row_title = "UP", row_title_rot = 0,
#                                 # column_order =order(colnames(df_enr_median_heat_UP)),
#                                 # row_order = order(rownames(df_enr_median_heat_UP)), 
#                                 show_row_names= TRUE,
#                                 show_column_names = FALSE,
#                                 # row_names_gp = grid::gpar(fontsize =8),
#                                 # column_names_gp = grid::gpar(fontsize =7),
#                                 # show_heatmap_legend = FALSE
#                                 top_annotation=ha,
#                                 # show_row_dend = TRUE,
#                                 # heatmap_legend_param = list(
#                                 # legend_direction = "horizontal", 
#                                 # legend_height = unit(1, "cm"),
#                                 # legend_gp = gpar(fontsize = 5)
#                                 )
# # combined gene and splice heatmaps
# png(file=paste0(fig_output_path,"GSVA_heatmap.png"),
#     width = 60,
#     height    = 50,
#     units     = "cm",
#     res       = 1200)

# draw(ht, legend_grouping = "original",heatmap_legend_side="bottom", annotation_legend_side="bottom")
# dev.off()

##########################
# PCA of GSVA res
##########################

print(dim(df_enr_all_sets))
prcomp.out <- prcomp(as.data.frame(t(df_enr_all_sets)),
                   center = TRUE,
                   scale. = TRUE)
  
# print(head(prcomp.out ))
print(dim(prcomp.out ))

# Merge PCA results with metadata
df_PCA <- data.frame(x = prcomp.out$x[,1],  y = prcomp.out$x[,2])
rownames(df_PCA) <- colnames(df_enr_all_sets)
pca.out.merge = cbind(df_PCA, df_meta)

print(dim(pca.out.merge ))

# Plot PCA all s
plt <- ggplot(pca.out.merge, aes(x, y, color = LUAD.Sample.Type)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="top",legend.title = element_blank()) +
    # scale_color_manual(values=pal) +
    labs(title= "", sep = ' ', x="PC1",y= "PC2")

# Save plot
ggsave(file.path(fig_output_path,
                paste("LUAD_GSVA_PCA.png", sep = '.')),
        device = "png",
        width = 5, height = 4,
        dpi = 300)


##########################
# PCA of GSVA res
##########################
umap.out <- umap(as.data.frame(prcomp.out$x), n_neighbors = 5, learning_rate = 0.5, 
                            init = "random", min_dist = 1, spread = 5)
umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
rownames(umap.out) <- colnames(df_enr_all_sets)

# Merge UMAP results with metadata
umap.out.merge = cbind(umap.out, df_meta)

# Plot UMAP
plt_umap <- ggplot(umap.out.merge, aes(x, y, color = LUAD.Sample.Type)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="top",legend.title = element_blank())+
    # guides(color=guide_legend(nrow=3,byrow=TRUE)) +
#   scale_color_manual(values=pal) +
    # scale_color_viridis_d() +
    # scale_color_viridis_d(option="C")+
    # scale_colour_brewer(palette = "Set2") +
        # scale_color_manual(values=c("#006400",
        # "#9A6324",
        #     "#ff0000",
        #     "#ffd700",
        #     "#00ff00",
        #     "#00ffff",
        #     "#ff00ff",
        #     "#ffb6c1"))+
    labs(title= "", sep = ' ', 
        x="UMAP1",
        y="UMAP2") 

# Save plot
ggsave(file.path(fig_output_path,
                paste0("LUAD_GSVA_UMAP.png")),
        device = "png",
        width = 5, height = 4,
        dpi = 300)
