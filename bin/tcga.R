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

# cat("are these weird...... ?")
# print(head(colnames(df_LUAD)))

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
df_LUAD_dcc <- read.table(file = "/mnt/tcga/dcc_sample_data_scrape_all.tsv",
                          sep="\t", header = TRUE) %>%
                    filter(cases.0.project.project_id == "TCGA-LUAD") 
df_tcga_clinical <- read.table(file = "/mnt/tcga/combined_study_clinical_data.tsv",
                          sep="\t", header = TRUE) 

df_meta <- merge(df_LUAD_dcc, df_tcga_clinical, by.x = "cases.0.submitter_id", by.y = "Patient.ID") 
# %>%
    # tibble::column_to_rownames("id")

cat("Dim of DCC meta")
print(dim(df_LUAD_dcc))

cat("Dim of clinical meta")
print(dim(df_tcga_clinical))

cat("Dim merged metadata ")
dim(df_meta)

intersect(colnames(df_LUAD_spliceset), df_meta$id)

print(head(colnames(df_LUAD_spliceset)))
print(head(df_meta$id))


print(colnames(df_LUAD_spliceset))
quit()



# add sample name to col names 
names(df_LUAD_meta) <- df_LUAD_meta[5, ]

# print(head(df_LUAD_meta))
print(dim(df_LUAD_meta))
print(row.names(df_LUAD_meta))

quit()

#########################
# Heatmap
#########################

df_LUAD_spliceset_clean <- df_LUAD_spliceset[which(rowMeans(!is.na(df_LUAD_spliceset)) > 0.5), ]
# print(dim(df_LUAD_spliceset_clean))
df_LUAD_spliceset_clean_scaled = t(scale(t(df_LUAD_spliceset_clean)))
print(head(df_LUAD_spliceset_clean_scaled))

# # calculate variance per column
# variances <- apply(df_LUAD_spliceset_clean[, -1], 1, var)
# print(head(variances))
# # sort variance, grab index of the first 2
# sorted <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:101] # replace 2 with 100 ...

# # use that to subset the original data
# df_LUAD_spliceset_highvar <- df_LUAD_spliceset_clean[, sorted]
# head(df_LUAD_spliceset_highvar)
# dim(df_LUAD_spliceset_highvar)

df_LUAD_spliceset_clean_sample <- df_LUAD_spliceset_clean_scaled[sample(nrow(df_LUAD_spliceset_clean_scaled), 20), ]
print(head(df_LUAD_spliceset_clean_sample))
print(dim(df_LUAD_spliceset_clean_sample))

# df_LUAD_spliceset_clean_sample <- t(df_LUAD_spliceset_clean_sample) %>% as.data.frame()

ht <- ComplexHeatmap::Heatmap(df_LUAD_spliceset_clean_sample,
                # cluster_rows = FALSE
                                # row_title = "UP", row_title_rot = 0,
                                # column_order =order(colnames(df_enr_median_heat_UP)),
                                # row_order = order(rownames(df_enr_median_heat_UP)), 
                                show_row_names= FALSE,
                                show_column_names = FALSE,
                                # row_names_gp = grid::gpar(fontsize =8),
                                # column_names_gp = grid::gpar(fontsize =7),
                                # show_heatmap_legend = FALSE
                                top_annotation=ha,
                                show_row_dend = FALSE,
                                # heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                # legend_height = unit(1, "cm"),
                                # legend_gp = gpar(fontsize = 5)
                                )
# combined gene and splice heatmaps
png(file=paste0(fig_output_path,"LUAD_spliceset_heatmap.png"),
    width = 40,
    height    = 20,
    units     = "cm",
    res       = 1200)

draw(ht)
dev.off()


#########################
# PCA
#########################