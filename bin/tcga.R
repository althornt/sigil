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

# df_LUAD_dcc <- read.table(file = "/mnt/tcga/dcc_sample_data_scrape_all.tsv",
#                           sep="\t", header = TRUE) %>%
#                     filter(cases.0.project.project_id == "TCGA-LUAD") 
df_tcga_clinical <- read.table(file = "/mnt/tcga/luad_tcga_pan_can_atlas_2018_clinical_data.tsv",
                          sep="\t", header = TRUE) %>%
                          filter(Patient.ID %in% df_bam_id$Case.ID)

df_meta <- merge(df_bam_id, df_tcga_clinical, by.x = "Case.ID", by.y = "Patient.ID", all= FALSE) 
# %>%
    # tibble::column_to_rownames("File.ID")

cat("Dim of bam meta")
print(dim(df_bam_id))

cat("Dim of clinical meta")
print(dim(df_tcga_clinical))

cat("Dim merged metadata ")
dim(df_meta)
# quit()

cat("Intersection of samples in metadata and PS table")
print(length(intersect(colnames(df_LUAD_spliceset), df_meta$File.Name)))
print(length(unique(df_meta$File.Name)))

print(colnames(df_meta))


df_meta <- df_meta %>%
    filter(File.Name %in% colnames(df_LUAD_spliceset)) %>%
    # distinct() %>%
    select(c("Tumor.Type", "File.Name", "MSI.MANTIS.Score", "Mutation.Count", "Person.Neoplasm.Cancer.Status", "Ragnum.Hypoxia.Score", "Sex", "Sample.Type.x")) %>%
    column_to_rownames("File.Name")
    # %>%
    # t() %>%
    # as.data.frame()

# head(df_meta)

# print(length(unique(df_meta["File.Name",])))

# names(df_meta) <- df_meta["File.Name",]
head(df_meta)
dim(df_meta)
# quit()

# ,"Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code"
# Tumor Type
# Tumor Stage 2009
# Primary Tumor
# Person Neoplasm Status
# Tissue Source Site	TMB (nonsynonymous)	Patient Smoking History Category
# Sex

#########################
# Heatmap
#########################

df_LUAD_spliceset_clean <- df_LUAD_spliceset[which(rowMeans(!is.na(df_LUAD_spliceset)) > 0.5), ] %>%
    select(rownames(df_meta))

# Drop rows with low variance
print(dim(df_LUAD_spliceset_clean))
getVar <- apply(df_LUAD_spliceset_clean[, -1], 1, var, na.rm=TRUE)
param <- 0.04
df_LUAD_spliceset_clean_highVar <- df_LUAD_spliceset_clean[getVar > param & !is.na(getVar), ]

print(dim(df_LUAD_spliceset_clean_highVar))

# Scale data 
df_LUAD_spliceset_clean_scaled = t(scale(t(df_LUAD_spliceset_clean_highVar)))
# print(head(df_LUAD_spliceset_clean_scaled))

ha <- HeatmapAnnotation(
    df = df_meta %>% select("Sample.Type.x"), 
    col = list(Sample.Type.x = c("Primary Tumor" ="green", "Solid Tissue Normal"= "orange", "Recurrent Tumor"="purple"    ))
)

ht <- ComplexHeatmap::Heatmap(df_LUAD_spliceset_clean_scaled,
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
    width = 50,
    height    = 20,
    units     = "cm",
    res       = 1200)

draw(ht)
dev.off()


#########################
# PCA
#########################


make_umap <- function(num_neighbor,meta_col,df_PCA,out_path) {

  set.seed(123)

  # Make color palette
  n <- length(unique(df_meta[[meta_col]]))
  print(unique(df_meta[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(df_PCA, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(df_PCA)

  # Merge UMAP results with metadata
  umap.out.merge = cbind(umap.out, df_meta)

  # Plot UMAP
  plt <- ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    guides(color = guide_legend(nrow = 4)) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP MESA: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(fig_output_path,
                   paste(out_path,meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)
}


# Only keep events over .75 variance cut off
var <- apply(df_LUAD_spliceset_clean[, -1], 1, var)
param <- quantile(var, c(.75), na.rm=T)
df_LUAD_spliceset_clean_filt <- df_LUAD_spliceset_clean[var > param & !is.na(var), ]

# Transpose and format
df_LUAD_spliceset_clean_filt_t <- as.data.frame(t(df_LUAD_spliceset_clean_filt))
rownames(df_LUAD_spliceset_clean_filt_t) <- colnames(df_LUAD_spliceset_clean)


# Make color palette
n <- length(unique(df_meta[["Sample.Type.x"]]))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# PCA
prcomp.out = as.data.frame(prcomp(na.omit(df_LUAD_spliceset_clean_filt_t), center=T,  scale = T)$x)

# Merge PCA results with metadata
df_PCA <- data.frame(x = prcomp.out[,1],  y = prcomp.out[,2])
pca.out.merge = cbind(df_PCA, df_meta)

# Plot PCA
plt <- ggplot(pca.out.merge, aes(x, y, color = Sample.Type.x)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= "", sep = ' ')

# Save plot
ggsave(file.path(fig_output_path,
                paste("PCA.png", sep = '.')),
        device = "png",
        width = 5, height = 4,
        dpi = 300)


# Plot UMAPS from PCA
# Making variations of UMAPs with different numbers of neighbors
# lapply(c(3,5,10, 20), make_umap, meta_col="Tumor.Type",
#     df_PCA = prcomp.out, out_path = paste0("PCA.UMAP"))
# lapply(c(3,5,10, 20), make_umap, meta_col="Person.Neoplasm.Cancer.Status",
#     df_PCA = prcomp.out, out_path = paste0("PCA.UMAP"))
# lapply(c(5, 20), make_umap, meta_col="Sex",
#     df_PCA = prcomp.out, out_path = paste0("PCA.UMAP"))
lapply(c(3,5,10, 20), make_umap, meta_col="Sample.Type.x",
    df_PCA = prcomp.out, out_path = paste0("PCA.UMAP"))


# Plot UMAPS from df
# Making variations of UMAPs with different numbers of neighbors
# lapply(c(3,5,10, 20), make_umap, meta_col="Tumor.Type",
#     df_PCA = df_LUAD_spliceset_clean_filt_t, out_path = paste0("UMAP"))
# lapply(c(3,5,10, 20), make_umap, meta_col="Person.Neoplasm.Cancer.Status",
#     df_PCA = df_LUAD_spliceset_clean_filt_t, out_path = paste0("UMAP"))
# lapply(c(5, 20), make_umap, meta_col="Sex",
#     df_PCA = df_LUAD_spliceset_clean_filt_t, out_path = paste0("UMAP"))
lapply(c(3,5,10, 20), make_umap, meta_col="Sample.Type.x",
    df_PCA = df_LUAD_spliceset_clean_filt_t, out_path = paste0("UMAP"))
