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
library(tidyr)

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
      "ls_mesa_cluster_files"=res_cluster_path
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

plot_PS_hist <- function(df, out_path){
  data_long <- df %>%
    tibble::rownames_to_column("event") %>%
    tidyr::gather(., key="sample", value = "PS", -c(event)) %>%
    as.data.frame()
  p <- ggplot( data_long, aes(x = PS)) +
      geom_histogram() +
      scale_y_continuous(trans='log2')
  ggsave(plot = p, filename = out_path)
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


# Open indiviidual count files
print(ls_mesa_inc_count_files)

# df_Song_allPS <- read.table("/mnt/results/sigil_results_SRP253519_Song_20211114/mesa_out/mesa_allPS.tsv",
#                               row.names = 1, header=T)
#
# print(dim(df_Song_allPS))
# print(max(df_Song_allPS))
# plot_PS_hist(df_Song_allPS,  paste0(opt$out_dir,"/df_Song_allPS.png"))

#
# print(dim(df_Song_inclusionCounts))
# print(max(df_Song_inclusionCounts))

# df_Song_inclusionCounts <- read.table("/mnt/results/sigil_results_SRP253519_Song_20211114/mesa_out/mesa_inclusionCounts.tsv",
#                               row.names = 1, header=T)
#
# print(dim(df_Song_inclusionCounts))
# print(max(df_Song_inclusionCounts))


# plot_PS_hist(df_Song_inclusionCounts,  paste0(opt$out_dir,"/df_Song_inclusionCounts.png"))

# df_Monaco_inclusionCounts <- read.table("/mnt/results/sigil_results_SRP125125_Monaco_20211109/mesa_out/mesa_inclusionCounts.tsv",
#                               row.names = 1, header=T)
#
# print(dim(df_Monaco_inclusionCounts))
# print(max(df_Monaco_inclusionCounts))

# plot_PS_hist(df_Monaco_inclusionCounts,  paste0(opt$out_dir,"/df_Monaco_inclusionCounts.png"))


# df_Choi_inclusionCounts <- read.table("/mnt/results/sigil_results_SRP150419_Choi_20211124/mesa_out/mesa_inclusionCounts.tsv",
#                               row.names = 1, header=T)
#
# print(dim(df_Choi_inclusionCounts))
# print(max(df_Choi_inclusionCounts))

# plot_PS_hist(df_Choi_inclusionCounts,  paste0(opt$out_dir,"/df_Choi_inclusionCounts.png"))

############################
# Merging mesa files
#############################
# Create manifest listing each inclusion count file
# write.table(as.data.frame(unlist(ls_mesa_inc_count_files)),
#             file=paste0(opt$out_dir,"/manifest_mesa_inclusionCounts_files.tsv"),
#             col.names = FALSE,
#             row.names = FALSE,
#             quote = FALSE)
# # Merge all inclusion count files siwth MESA select_samples command
# # 2>&1 sends standard error standard output
# inc_cmd <- paste0(
#     "mesa select -m ",
#     opt$out_dir,"/manifest_mesa_inclusionCounts_files.tsv -o ",
#     opt$out_dir,"/merged_mesa_inclusionCounts.tsv  --join intersection 2>&1"
#     )
#
# print(inc_cmd)
# system(inc_cmd)

# Create manifest listing each PS file
# write.table(as.data.frame(unlist(ls_mesa_allPS_files)),
#             file=paste0(opt$out_dir,"/manifest_mesa_allPS_files.tsv"),
#             col.names = FALSE,
#             row.names = FALSE,
#             quote = FALSE)
#
# # Merge all PS file siwth MESA select_samples command
# # 2>&1 sends standard error standard output
# PS_cmd <- paste0(
#     "mesa select -m ",
#     opt$out_dir,"/manifest_mesa_allPS_files.tsv -o ",
#     opt$out_dir,"/merged_mesa_allPS.tsv  --join intersection 2>&1"
#     )
#
# print(PS_cmd)
# system(PS_cmd)


######################################
# UMAPs of PS before batch correction
######################################
# Read in merged allPS file
df_merged_allPS <- read.table(paste0(opt$out_dir,"/merged_mesa_allPS.tsv"),
                              row.names = 1, header=T)
# print(head(df_merged_allPS))
print(dim(df_merged_allPS))
plot_PS_hist(df_merged_allPS,  paste0(opt$out_dir,"/hist_allPS.png"))
print("sum")
print(sum(rowSums(is.na(df_merged_allPS))))

# Drop genes with low variance.
allPS_var <- apply(df_merged_allPS[, -1], 1, var)
print(length(allPS_var))

# For gene I used median (50% quantile) as the cutoff
# For splicing using 75% due to there being many more rows)
PS_param <- quantile(allPS_var, c(.75), na.rm=T)
print(PS_param)

df_allPS_filt <- df_merged_allPS[allPS_var > PS_param & !is.na(allPS_var), ]
print(dim(df_allPS_filt))

# Transpose and format
df_allPS_filt_t <- as.data.frame(t(df_allPS_filt))
rownames(df_allPS_filt_t) <- colnames(df_allPS_filt)

print(dim(df_allPS_filt))
print(dim(df_allPS_filt_t))
# PCA.
all.ps.prcomp.out = as.data.frame(prcomp(na.omit(df_allPS_filt_t), center=T,  scale = T)$x)

# Making variations of UMAPs with different numbers of neighbors
lapply(c(20), make_umap, meta_col="data_source",
  df_PCA = all.ps.prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")
lapply(c(20), make_umap, meta_col="LM22",
  df_PCA = all.ps.prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")
lapply(c(20), make_umap, meta_col="sigil_general",
  df_PCA = all.ps.prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")
lapply(c(20), make_umap, meta_col="LM6",
  df_PCA = all.ps.prcomp.out, out_path = "UMAPs_pre_batch_correction/mesa_incl_count_PCA_UMAP")


########################################################
# Batch correction of inclusion counts
#########################################################

# Read in merged inclusion count file
df_merged_inc_counts <- read.table(paste0(opt$out_dir,"/merged_mesa_inclusionCounts.tsv"),
                                    row.names = 1, header=T)

print("df_merged_inc_counts")
print(sum((is.na(df_merged_inc_counts))))
plot_PS_hist(df_merged_inc_counts,  paste0(opt$out_dir,"/merged_inc_counts.png"))

# Log2 + 1 transform counts
df_log2trans_inc_counts <- as.data.frame(log2(df_merged_inc_counts +1))
print(dim(df_log2trans_inc_counts))

plot_PS_hist(df_log2trans_inc_counts,  paste0(opt$out_dir,"/hist_inc_counts_log.png"))

print("df_log2trans_inc_counts")
print(sum((is.na(df_log2trans_inc_counts))))

df_mesa_inc_count_merge_log2_batch_corr <- limma::removeBatchEffect(
                                  df_log2trans_inc_counts,
                                  batch = df_merged_metadata$data_source,
                                  batch2 = df_merged_metadata$type
                                  )

print("df_mesa_inc_count_merge_log2_batch_corr")
print(dim(df_mesa_inc_count_merge_log2_batch_corr))
print(sum((is.na(df_mesa_inc_count_merge_log2_batch_corr))))

# Undo log2(x+1) with 2^x - 1
df_mesa_inc_count_merge_bc_counts = 2^df_mesa_inc_count_merge_log2_batch_corr -1
rownames(df_mesa_inc_count_merge_bc_counts) <- rownames(df_merged_inc_counts)
print(typeof(df_mesa_inc_count_merge_bc_counts))

plot_PS_hist(as.data.frame(df_mesa_inc_count_merge_bc_counts),
              paste0(opt$out_dir,"/hist_inc_counts_bc.png"))

# print(head(df_merged_inc_counts))
# print(head(df_mesa_inc_count_merge_bc_counts))


# df_mesa_inc_count_merge_bc_counts_r <- round(df_mesa_inc_count_merge_bc_counts, digits=0)

# print(head(df_mesa_inc_count_merge_bc_counts))
# print(head(df_mesa_inc_count_merge_bc_counts_r))
#
# print(dim(df_mesa_inc_count_merge_bc_counts))
# print(dim(df_mesa_inc_count_merge_bc_counts_r))

write.table(
  df_mesa_inc_count_merge_bc_counts,
  file.path(opt$out_dir,"batch_corr_mesa_inclusionCounts.tsv"),
  sep="\t",quote=F, col.names = NA)


########################################################
# Convert batch corrected counts to PS
#########################################################
# Merge all PS files with MESA select_samples command
# 2>&1 sends standard error standard output
PS_cmd <- paste0(
    "mesa counts_to_ps -i ",
    opt$out_dir,"/batch_corr_mesa_inclusionCounts.tsv -o ",
    opt$out_dir,"/batch_corr_mesa_allPS.tsv  --recluster 2>&1"
    )
print(PS_cmd)
system(PS_cmd)


########################################################
# UMAPs of PS after batch correction
#########################################################
df_merged_allPS_bc <- read.table(paste0(opt$out_dir,"/batch_corr_mesa_allPS.tsv"),
                              row.names = 1, header=T)
print("df_merged_allPS_bc")
print(head(df_merged_allPS_bc))
print(dim(df_merged_allPS_bc))
plot_PS_hist(df_merged_allPS_bc,  paste0(opt$out_dir,"/hist_allPS_bc.png"))

# Log2 + 1 transform PS
df_allPS_log2trans_bc <- as.data.frame(log2(df_merged_allPS_bc +1))
print("df_allPS_log2trans_bc")
print(dim(df_allPS_log2trans_bc))
plot_PS_hist(df_allPS_log2trans_bc,  paste0(opt$out_dir,"/hist_allPS_bc_log.png"))

print("-------------")
# # Drop genes with low variance.
allPS_log2trans_var_bc <- apply(df_allPS_log2trans_bc[, -1], 1, var)
print(length(allPS_log2trans_var_bc))

# For gene I used median (50% quantile) as the cutoff
# For splicing using 75% due to there being many more rows)
PS_param_bc <- quantile(allPS_log2trans_var_bc, c(.75), na.rm=T)
print(PS_param_bc)

df_allPS_log2trans_bc_filt <- df_allPS_log2trans_bc[allPS_log2trans_var_bc > PS_param_bc & !is.na(allPS_log2trans_var_bc), ]
print(dim(df_allPS_log2trans_bc_filt))

# # Transpose and format
df_allPS_log2trans_bc_filt_t <- as.data.frame(t(df_allPS_log2trans_bc_filt))
rownames(df_allPS_log2trans_bc_filt_t) <- colnames(df_allPS_log2trans_bc_filt)

print(dim(df_allPS_log2trans_bc_filt))
print(dim(df_allPS_log2trans_bc_filt_t))

# PCA.
all.ps.prcomp.out.bc = as.data.frame(prcomp(na.omit(df_allPS_log2trans_bc_filt_t), center=T,  scale = T)$x)

# Making variations of UMAPs with different numbers of neighbors
lapply(c(10,20,30), make_umap, meta_col="data_source",
  df_PCA = all.ps.prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
lapply(c(10,20,30), make_umap, meta_col="LM22",
  df_PCA = all.ps.prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
lapply(c(10,20,30), make_umap, meta_col="sigil_general",
  df_PCA = all.ps.prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")
lapply(c(10,20,30), make_umap, meta_col="LM6",
  df_PCA = all.ps.prcomp.out.bc, out_path = "UMAPs_post_batch_correction/mesa_incl_count_PCA_UMAP")

# # Drop non LM22 samples from mesa PS
df_merged_allPS_bc_lm22 <- df_merged_allPS_bc %>%
  dplyr::select(ls_smpls_lm22)
rownames(df_merged_allPS_bc_lm22) <- rownames(df_merged_allPS_bc)
print(head(df_merged_allPS_bc_lm22))
write.table(
  df_merged_allPS_bc_lm22,
  file.path(opt$out_dir,"batch_corr_mesa_allPS_LM22.tsv"), quote=F,sep="\t", na="nan",
  col.names = NA, row.names= TRUE)

print("df_merged_allPS_bc_lm22")
print(head(rownames(df_merged_allPS_bc_lm22)))
print(dim(df_merged_allPS_bc_lm22))



#################################
# Compare PS before and after
################################
