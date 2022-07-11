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

unique(colnames(df_splice_set_byevent))


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
    select(c("Tumor.Type", "File.Name", "MSI.MANTIS.Score", "Mutation.Count", 
                "Person.Neoplasm.Cancer.Status", "Ragnum.Hypoxia.Score", 
                "Sex", "Sample.Type.x", "Cancer.Type", "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code")) %>%
    column_to_rownames("File.Name") %>%
    rename(LUAD.Sample.Type = Sample.Type.x) %>%
    arrange(rownames(.))


head(df_meta)
# quit()
dim(df_meta)
df_bam_id <- df_bam_id %>%
    filter(File.Name %in% colnames(df_LUAD_spliceset)) %>%
    column_to_rownames("File.Name") %>%
    arrange(rownames(.))
head(df_bam_id)

head(df_meta)
# quit()
# ,"Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code"
# Tumor Type
# Tumor Stage 2009
# Primary Tumor
# Person Neoplasm Status
# Tissue Source Site	TMB (nonsynonymous)	Patient Smoking History Category
# Sex


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
# Counts of sets in the sigil LUAD junctions 
################################################

head(df_splice_set)

df_set_cnt <- df_splice_set %>%
    filter(event %in% rownames(df_mesa_tumor_vs_norm_sigil)) %>%
    count(set) %>%
    arrange(desc(n))

print(head(df_set_cnt, n=20))
print(dim(df_set_cnt))

# print(head(table(df_set_cnt)))


# ggplot(as.data.frame(tbl), aes(factor(Depth), Freq, fill = Species)) +     
#   geom_col(position = 'dodge')

#     ggplot(yes, aes(x=key, y=value)) + 
#   geom_bar(stat="identity")

################################################
# Heatmap sigil set
################################################
df_LUAD_spliceset_high_var <- df_LUAD_spliceset %>%
    select( intersect(colnames(.),rownames(df_meta))) %>%
    # filter(rownames(.) %in% rownames(df_mesa_tumor_vs_norm_sigil))%>%
    rowwise() %>%
    mutate(variance = c_across(everything()) %>% var()) 

rownames(df_LUAD_spliceset_high_var) <- rownames(df_LUAD_spliceset)
df_LUAD_spliceset_high_var <- df_LUAD_spliceset_high_var %>%
     ungroup() %>%
     arrange(desc(variance)) %>%
     head(300) %>%
    as.data.frame() %>%
    select(!variance) %>%
    select(sort(names(.)))

# print(head(df_LUAD_spliceset_high_var))
print(dim(df_LUAD_spliceset_high_var))
# quit()
df_LUAD_spliceset_scaled = t(scale(t(df_LUAD_spliceset_high_var)))

ha <- HeatmapAnnotation(
    df = df_meta %>% select(c("LUAD.Sample.Type")) , 
    col = list(LUAD.Sample.Type = c("Primary Tumor" ="#F8766D", 
                            "Solid Tissue Normal"= "#619CFF",
                             "Recurrent Tumor"="#00BA38"    )
                # Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code = c(
                #                     "NA" = "gray" ,
                #                     "STAGE IB" = "black" ,
                #                     "STAGE I" = "black" ,
                #                     "STAGE IA" = "black" ,
                #                     "STAGE II" = "yellow" ,
                #                     "STAGE IIA" = "yellow" ,
                #                     "STAGE IIB" = "yellow" ,
                #                     "STAGE IIIA" = "orange" ,
                #                     "STAGE IIIB" = "orange" ,
                #                     "STAGE IV" = "red" 
            
                )            )

ht <- ComplexHeatmap::Heatmap(df_LUAD_spliceset_scaled,
                    name = "Z-Score PS",
                # cluster_columns = FALSE,
                                # row_title = "UP", row_title_rot = 0,
                                # column_order =order(colnames(df_enr_median_heat_UP)),
                                # row_order = order(rownames(df_enr_median_heat_UP)), 
                                show_row_names= FALSE,
                                show_column_names = FALSE,
                                # row_names_gp = grid::gpar(fontsize =8),
                                # column_names_gp = grid::gpar(fontsize =7),
                                # show_heatmap_legend = FALSE
                                top_annotation=ha,
                                # show_row_dend = TRUE,
                                # heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                # legend_height = unit(1, "cm"),
                                # legend_gp = gpar(fontsize = 5)
                                )
# combined gene and splice heatmaps
png(file=paste0(fig_output_path,"LUAD_spliceset_heatmap.png"),
    width = 50,
    height    = 20,
    units     = "cm",
    res       = 1200)

draw(ht, legend_grouping = "original",heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()

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
                   scale. = TRUE)$x
  
# print(head(prcomp.out ))
print(dim(prcomp.out ))

# Merge PCA results with metadata
df_PCA <- data.frame(x = prcomp.out[,1],  y = prcomp.out[,2])
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
    labs(title= "", sep = ' ')

# Save plot
ggsave(file.path(fig_output_path,
                paste("LUAD_spliceset_PCA.png", sep = '.')),
        device = "png",
        width = 5, height = 4,
        dpi = 300)

############################################################################
################################################
# Heatmap tumor vs normal + sigil set + matched delta
################################################
df_LUAD_sig_diff_and_sigil <- df_LUAD_spliceset %>%
    select( intersect(colnames(.),rownames(df_meta))) %>%
    # filter(rownames(.) %in% rownames(df_mesa_tumor_vs_norm_sigil))%>%
    filter(rownames(.) %in% matches$event)%>%
    select(sort(names(.))) 

df_LUAD_sig_diff_and_sigil_scaled = t(scale(t(df_LUAD_sig_diff_and_sigil)))

ha <- HeatmapAnnotation(
    df = df_meta %>% select(c("LUAD.Sample.Type")), 
    col = list(LUAD.Sample.Type = c("Primary Tumor" ="#F8766D", 
                            "Solid Tissue Normal"= "#619CFF",
                             "Recurrent Tumor"="#00BA38"    )
                # Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code = c(
                #                     "NA" = "gray" ,
                #                     "STAGE IB" = "black" ,
                #                     "STAGE I" = "black" ,
                #                     "STAGE IA" = "black" ,
                #                     "STAGE II" = "yellow" ,
                #                     "STAGE IIA" = "yellow" ,
                #                     "STAGE IIB" = "yellow" ,
                #                     "STAGE IIIA" = "orange" ,
                #                     "STAGE IIIB" = "orange" ,
                #                     "STAGE IV" = "red" 
            
                )            )

ht <- ComplexHeatmap::Heatmap(df_LUAD_sig_diff_and_sigil_scaled,
                    name = "Z-Score PS",
                # cluster_columns = FALSE,
                                # row_title = "UP", row_title_rot = 0,
                                # column_order =order(colnames(df_enr_median_heat_UP)),
                                # row_order = order(rownames(df_enr_median_heat_UP)), 
                                show_row_names= FALSE,
                                show_column_names = FALSE,
                                # row_names_gp = grid::gpar(fontsize =8),
                                # column_names_gp = grid::gpar(fontsize =7),
                                # show_heatmap_legend = FALSE
                                top_annotation=ha,
                                # show_row_dend = TRUE,
                                heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                legend_height = unit(10, "cm"),
                                legend_gp = gpar(fontsize = 20)
                                ))
# combined gene and splice heatmaps
png(file=paste0(fig_output_path,"LUAD_sig_diff_sigil_heatmap_scaled.png"),
    width = 60,
    height    = 20,
    units     = "cm",
    res       = 1200)

draw(ht, legend_grouping = "original",heatmap_legend_side="right", annotation_legend_side="right")

ht <- ComplexHeatmap::Heatmap(df_LUAD_sig_diff_and_sigil,
                    name = "PS",
                      col = circlize::colorRamp2(c(0, 1.0), c("white", "red")),

                # cluster_columns = FALSE,
                                # row_title = "UP", row_title_rot = 0,
                                # column_order =order(colnames(df_enr_median_heat_UP)),
                                # row_order = order(rownames(df_enr_median_heat_UP)), 
                                show_row_names= FALSE,
                                show_column_names = FALSE,
                                # row_names_gp = grid::gpar(fontsize =8),
                                # column_names_gp = grid::gpar(fontsize =7),
                                # show_heatmap_legend = FALSE
                                top_annotation=ha,
                                # show_row_dend = TRUE,
                                heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                legend_height = unit(2, "cm"),
                                # legend_gp = gpar(fontsize = 5), 
                                legend_gp = gpar(fontsize = 14)
                                ))
# combined gene and splice heatmaps
png(file=paste0(fig_output_path,"LUAD_sig_diff_sigil_heatmap_unscaled.png"),
    width = 60,
    height    = 20,
    units     = "cm",
    res       = 1200)

draw(ht, legend_grouping = "original",heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

###################################################
# PCA tumor vs normal + sigil set
####################################################

df_LUAD_sig_diff_and_sigil_clean <- df_LUAD_sig_diff_and_sigil[which(rowMeans(!is.na(df_LUAD_sig_diff_and_sigil)) > 0.5), ]  %>%
    select(intersect(colnames(df_LUAD_sig_diff_and_sigil),rownames(df_meta))) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = T), x))

print(dim(df_LUAD_sig_diff_and_sigil))

prcomp.out <- prcomp(as.data.frame(t(df_LUAD_sig_diff_and_sigil_clean)),
                   center = TRUE,
                   scale. = TRUE)$x
  
# print(head(prcomp.out ))
print(dim(prcomp.out ))

# Merge PCA results with metadata
df_PCA <- data.frame(x = prcomp.out[,1],  y = prcomp.out[,2])
rownames(df_PCA) <- colnames(df_LUAD_sig_diff_and_sigil_clean)
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
    labs(title= "", sep = ' ')

# Save plot
ggsave(file.path(fig_output_path,
                paste("LUAD_sig_diff_sigil_PCA.png", sep = '.')),
        device = "png",
        width = 5, height = 4,
        dpi = 300)

