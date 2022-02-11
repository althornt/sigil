#!/usr/bin/env Rscript

library(optparse)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(uwot)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(ensembldb)
library(DESeq2)
library(pheatmap)
library(limma)
library(purrr)

importMetaMESA <- function(row){
  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general, LM22, LM6)

  # Add metadata to column
  df_metadata$data_source <- row[1] # add name of data source
  df_metadata$type <- row[4] # add rna-seq type (paired vs single)

  # Get paths to MESA inclusion count files and allPS files
  res_inc_count_path <- file.path(row[2], "mesa_out", "mesa_inclusionCounts.tsv")
  res_allPS_path <- file.path(row[2], "mesa_out", "mesa_allPS.tsv")
  res_cluster_path <- file.path(row[2], "mesa_out", "mesa_allClusters.tsv")

  return(list(
      "ls_mesa_inc_count_files"=res_inc_count_path,
      "ls_mesa_allPS_files"=res_allPS_path,
      "ls_mesa_cluster_files"=res_cluster_path,
      "metadata"=df_metadata,
      "ls_samples_run"=df_metadata$Run))
}

make_umap <- function(num_neighbor,meta_col,df_PCA,out_path) {

  set.seed(123)

  # Make color palette
  n <- length(unique(df_merged_metadata_lm22[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(df_PCA, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(df_PCA)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, df_merged_metadata_lm22)

  # Plot UMAP
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP MESA: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste(out_path,meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)
}

###################
# MAIN
###################
# Arguments
option_list <- list(

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "path to write outputs"),

  optparse::make_option(
    c("-m", "--manifest"),
    type = "character",
    default = NULL,
    help = "path to tsv manifest file with 3 columns: data_source, res_path, metadata_path, type")

  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/UMAPs_pre_batch_correction/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs_pre_batch_correction/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/UMAPs_post_batch_correction/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs_post_batch_correction/"),
              recursive = TRUE, showWarnings = TRUE)}

# Open manifest
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(df_manifest)

# Import and combine source metadata files
ls_mesa_meta = apply(df_manifest, 1, importMetaMESA)

# Split into mesa and metadata files for each data set
ls_mesa_inc_count_files <- ls_mesa_allPS_files <- ls_mesa_cluster_files <- ls_meta <- ls_sample_names <- c()
for (item in ls_mesa_meta) {
     ls_mesa_inc_count_files <- append(ls_mesa_inc_count_files, item[1])
     ls_mesa_allPS_files <- append(ls_mesa_allPS_files, item[2])
     ls_mesa_cluster_files <- append(ls_mesa_cluster_files, item[3])
     ls_meta <- append(ls_meta, item[4])
     ls_sample_names <- append(ls_sample_names, item[5])
   }

# Combine metadata from each data source by rows
df_merged_metadata <- do.call("rbind", ls_meta)
rownames(df_merged_metadata) <- c()

# Remove samples without LM22 labels from metadata
df_merged_metadata_lm22 <- df_merged_metadata %>%
   dplyr::filter(LM22 != "")
write.csv(df_merged_metadata_lm22,
            file.path(file.path(opt$out_dir,"lm22_metadata.csv")),
            row.names = FALSE, quote=F)

# List of samples with LM22 label
ls_smpls_lm22 <- as.character(df_merged_metadata_lm22$Run)

############################
# Merging mesa files
#############################
# Create manifest listing each inclusion count file
write.table(as.data.frame(unlist(ls_mesa_inc_count_files)),
            file=paste0(opt$out_dir,"/manifest_mesa_inclusionCounts_files.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
# Merge all inclusion count files siwth MESA select_samples command
# 2>&1 sends standard error standard output
inc_cmd <- paste0(
    "mesa select_samples -m ",
    opt$out_dir,"/manifest_mesa_inclusionCounts_files.tsv -o ",
    opt$out_dir,"/merged_mesa_inclusionCounts.tsv  --join intersection 2>&1"
    )

print(inc_cmd)
system(inc_cmd)

# Create manifest listing each PS file
write.table(as.data.frame(unlist(ls_mesa_allPS_files)),
            file=paste0(opt$out_dir,"/manifest_mesa_allPS_files.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

# Merge all PS file siwth MESA select_samples command
# 2>&1 sends standard error standard output
PS_cmd <- paste0(
    "mesa select -m ",
    opt$out_dir,"/manifest_mesa_allPS_files.tsv -o ",
    opt$out_dir,"/merged_mesa_allPS.tsv  --join intersection 2>&1"
    )

print(PS_cmd)
system(PS_cmd)


######################################
# UMAPs of PS before batch correction
######################################

# Read in merged allPS file
# df_merged_allPS <- read.table(paste0(opt$out_dir,"/merged_mesa_allPS.tsv"),
#                               row.names = 1, header=T)
# print(head(df_merged_allPS))
# print(head(rownames(df_merged_allPS)))
#
# # Log2 + 1 transform PS
# log2trans_dat <- as.data.frame(log2(df_merged_allPS +1))
# print(dim(log2trans_dat))


# # # Drop genes with low variance.
# # getVar <- apply(log2trans_dat[, -1], 1, var)
# print(dim(getVar))








# # For gene I used median (50% quantile) as the cutoff
# # For splicing using 75% due to there being many more rows)
# param <- quantile(getVar, c(.75))
# log2trans_dat_filt <- log2trans_dat[getVar > param & !is.na(getVar), ]
#
# # Transpose and format
# log2trans_dat_filt_t <- as.data.frame(t(log2trans_dat_filt))
# rownames(log2trans_dat_filt_t) <- colnames(log2trans_dat_filt)
#
# # PCA.
# prcomp.out = as.data.frame(prcomp(log2trans_dat_filt_t, scale = F)$x)
#
# # Making variations of UMAPs with different numbers of neighbors
# lapply(c(20), make_umap, meta_col="data_source",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")
# lapply(c(20), make_umap, meta_col="LM22",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")
# lapply(c(20), make_umap, meta_col="sigil_general",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")
# lapply(c(20), make_umap, meta_col="LM6",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")


########################################################
# Batch correction of inclusion counts
#########################################################

# Read in merged inclusion count file
# df_merged_inc_counts <- read.table(paste0(opt$out_dir,"/merged_mesa_inclusionCounts.tsv"))

# Log2 + 1 transform counts
# log2trans_dat <- as.data.frame(log2(df_merged_allPS +1))
# print(dim(log2trans_dat))


# df_mesa_inc_count_merge_log2_batch_corr <- limma::removeBatchEffect(
#                                   log2trans_dat,
#                                   batch = df_merged_metadata$data_source,
#                                   batch2 = df_merged_metadata$type
#                                   )

########################################################
# Convert batch corrected counts to PS
#########################################################



########################################################
# UMAPs of PS after batch correction
#########################################################



# #######################
# # Batch correction
# #######################
# df_mesa_inc_count_merge_log2_batch_corr <- limma::removeBatchEffect(
#                                   log2trans_dat,
#                                   batch = df_merged_metadata$data_source,
#                                   batch2 = df_merged_metadata$type
#                                   )
#
# # Drop genes with low variance.
# getVar_bc <- apply(df_mesa_inc_count_merge_log2_batch_corr[, -1], 1, var)
#
# # For gene I used median (50% quantile) as the cutoff
# # For splicing using 75% due to therembeing much more rows)
# param_bc <- quantile(getVar_bc, c(.75))
# log2trans_dat_filt_bc <- df_mesa_inc_count_merge_log2_batch_corr[getVar_bc > param_bc & !is.na(getVar_bc), ]
#
# Transpose and format
# log2trans_dat_filt_t_bc <- as.data.frame(t(log2trans_dat_filt_bc))
# rownames(log2trans_dat_filt_t_bc) <- colnames(log2trans_dat_filt_bc)
#
# # PCA.
# prcomp.out.bc = as.data.frame(prcomp(log2trans_dat_filt_t_bc, scale = F)$x)
#
# # Making variations of UMAPs with different numbers of neighbors
# lapply(c(20), make_umap, meta_col="data_source",
#   df_PCA = prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
# lapply(c(20), make_umap, meta_col="LM22",
#   df_PCA = prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
# lapply(c(20), make_umap, meta_col="sigil_general",
#   df_PCA = prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
# lapply(c(20), make_umap, meta_col="LM6",
#   df_PCA = prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
##############################################
# Undo log2(x+1) with 2^x - 1
##############################################
# df_mesa_inc_count_merge_bc_counts = 2^df_mesa_inc_count_merge_log2_batch_corr -1
# rownames(df_mesa_inc_count_merge_bc_counts) <- rownames(df_mesa_inc_count_merge)
#
# # print(head(df_mesa_inc_count_merge))
# # print(head(df_mesa_inc_count_merge_bc_counts))
#
#
# df_mesa_inc_count_merge_bc_counts_r <- round(df_mesa_inc_count_merge_bc_counts, digits=0)
#
# # print(head(df_mesa_inc_count_merge_bc_counts))
# # print(head(df_mesa_inc_count_merge_bc_counts_r))
# #
# # print(dim(df_mesa_inc_count_merge_bc_counts))
# # print(dim(df_mesa_inc_count_merge_bc_counts_r))
#
# write.table(
#   df_mesa_inc_count_merge_bc_counts_r,
#   file.path(opt$out_dir,"batch_corr_mesa_inclusionCounts.tsv"),
#   sep="\t",quote=F, col.names = NA)
# ##############################################
# # MESA quant on batch corrected counts
# ##############################################
#
# # MESA counts_to_ps command, 2>&1 sends standard error standard output
# cmd <- paste0(
#     "mesa counts_to_ps -i ",
#     opt$out_dir,"/batch_corr_mesa_inclusionCounts.tsv -c ",
#     opt$out_dir,"/mesa_allClusters.tsv -o ",
#     opt$out_dir,"/batch_corr_mesa_allPS.tsv  2>&1"
# )
#
# print(cmd)
# system(cmd)
#
#
#
# df_bc_allPS <- read.table(paste0(opt$out_dir,"/batch_corr_mesa_allPS.tsv"))
# print(dim(df_bc_allPS))

#
# # Remove samples without LM22 labels from metadata
# df_merged_metadata_lm22 <- df_merged_metadata %>%
#    dplyr::filter(LM22 != "")
# write.csv(df_merged_metadata_lm22,
#             file.path(file.path(opt$out_dir,"lm22_metadata.csv")),
#             row.names = FALSE, quote=F)
#
# # List of samples with LM22 lab
# ls_smpls_lm22 <- as.character(df_merged_metadata_lm22$Run)
#
# # Drop non LM22 samples from mesa counts
# df_mesa_inc_count_merge_lm22 <- df_mesa_inc_count_merge %>%
#   dplyr::select(ls_smpls_lm22)
# rownames(df_mesa_inc_count_merge_lm22) <- df_mesa_inc_count_merge$cluster
# write.table(
#   df_mesa_inc_count_merge_lm22,
#   file.path(opt$out_dir,"LM22_mesa_inclusionCounts.tsv"),quote=F,sep="\t", col.names = NA)
# #
# # print("df_mesa_inc_count_merge_lm22")
# # print(dim(df_mesa_inc_count_merge_lm22))
#
#
# # # Drop non LM22 samples from mesa PS
# df_mesa_allPS_merge_lm22 <- df_mesa_allPS_merge %>%
#   dplyr::select(ls_smpls_lm22)
# rownames(df_mesa_allPS_merge_lm22) <- df_mesa_allPS_merge$cluster
# write.table(
#   df_mesa_allPS_merge_lm22,
#   file.path(opt$out_dir,"LM22_mesa_allPS.tsv"), quote=F,sep="\t", na="nan", col.names = NA)
#
# # print("df_mesa_allPS_merge_lm22")
# # print(dim(df_mesa_allPS_merge_lm22))



######### old


# # Merge mesa PS files by common clusters
# df_mesa_allPS_merge <- unlist(ls_mesa_allPS_files) %>%
#   lapply(read.csv, sep = "\t") %>%
#   purrr::reduce(inner_join, by = "cluster")
# write.table(
#   df_mesa_allPS_merge,
#   file.path(opt$out_dir,"mesa_allPS.tsv"), quote=F,sep="\t", na="nan", col.names = NA)
# print("df_mesa_allPS_merge")
# print(dim(df_mesa_allPS_merge))
#
# # Merge mesa inclusion count files
# df_mesa_inc_count_merge <- unlist(ls_mesa_inc_count_files) %>%
#   lapply(read.csv, sep = "\t") %>%
#   purrr::reduce(inner_join, by = "cluster")
# write.table(
#   df_mesa_inc_count_merge,
#   file.path(opt$out_dir,"mesa_inclusionCounts.tsv"),quote=F, sep="\t", col.names = NA)
# print("df_mesa_inc_count_merge")
# print(dim(df_mesa_inc_count_merge))


# Merge mesa cluster files and remove duplicate clusters
# AND filter out of first fol and secold column which is a list
# ls_df_mesa_clusters <- unlist(ls_mesa_cluster_files) %>%
#   lapply(read.table, sep = "\t", header= F)
#
# df_mesa_clusters_merge <- do.call("rbind", ls_df_mesa_clusters) %>%
#   distinct(.)
#
#
# print(head(df_mesa_clusters_merge))
# print(dim(df_mesa_clusters_merge))
#
# df_mesa_clusters_merge <- df_mesa_clusters_merge %>%
#   dplyr::filter(V1 %in% df_mesa_inc_count_merge$cluster)
#
# print(head(df_mesa_clusters_merge))
#
# df_mesa_clusters_merge <- df_mesa_clusters_merge %>%
#   dplyr::mutate(V2 = map(V2, ~ purrr::keep(.x, .x %in% df_mesa_inc_count_merge$cluster)))
#
# print(head(df_mesa_clusters_merge))
# print(dim(df_mesa_clusters_merge))
#
# print("-------")
#
# print(head(df_mesa_inc_count_merge$cluster))
#
# write.table(
#   df_mesa_clusters_merge,
#   file.path(opt$out_dir,"mesa_allClusters.tsv"),
#   quote=F,sep="\t", col.names = FALSE, row.names = FALSE)
#
# print("df_mesa_clusters_merge")
# print(dim(df_mesa_clusters_merge))


# Format merged mesa inclusion count files for log2 transformation
# rownames(df_mesa_inc_count_merge) <- df_mesa_inc_count_merge$cluster
# df_mesa_inc_count_merge <- df_mesa_inc_count_merge %>%
#   dplyr::select( -one_of(c("cluster")))
