#!/usr/bin/env Rscript

# Combine mesa files from different mesa runs and perform batch correction

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
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)


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
  res_IR_table_path <- file.path(row[2], "mesa_out", "mesa_ir_table_intron_retention.tsv")
  res_IR_cov_dir_path <- file.path(row[2], "mesa_out", "mesa_intron_coverage")


  return(list(
      "ls_mesa_inc_count_files"=res_inc_count_path,
      "ls_mesa_allPS_files"=res_allPS_path,
      "ls_mesa_cluster_files"=res_cluster_path,
      "ls_mesa_IR_table_files" = res_IR_table_path,
      "ls_mesa_IR_cov_dir" = res_IR_cov_dir_path,
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
  plt <- ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
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

df_to_UMAP <- function(input_df, output_path){

  # Only keep events over .75 variance cut off
  var <- apply(input_df[, -1], 1, var)
  param <- quantile(var, c(.75), na.rm=T)
  input_df_filt <- input_df[var > param & !is.na(var), ]

  # Transpose and format
  input_df_filt_t <- as.data.frame(t(input_df_filt))
  rownames(input_df_filt_t) <- colnames(input_df)

  # PCA
  prcomp.out = as.data.frame(prcomp(na.omit(input_df_filt_t), center=T,  scale = T)$x)

  # Making variations of UMAPs with different numbers of neighbors
  lapply(c(20), make_umap, meta_col="data_source",
    df_PCA = prcomp.out, out_path = paste0(output_path))
  lapply(c(20), make_umap, meta_col="LM22",
    df_PCA = prcomp.out, out_path = paste0(output_path))
  lapply(c(20), make_umap, meta_col="sigil_general",
    df_PCA = prcomp.out, out_path = paste0(output_path))
  lapply(c(20), make_umap, meta_col="LM6",
    df_PCA = prcomp.out, out_path = paste0(output_path))
}

plot_before_after <- function(df_before, df_after, filename, val){

  before_long <- df_before %>%
    tibble::rownames_to_column("event") %>%
    tidyr::gather(., key="sample", value = val, -c(event)) %>%
    as.data.frame()

  after_long <- df_after %>%
    tibble::rownames_to_column("event") %>%
    tidyr::gather(., key="sample", value = val, -c(event)) %>%
    as.data.frame()

  p <- ggplot() +
      geom_freqpoly(data = before_long, aes(x = val), colour =  "orange") +
      geom_freqpoly(data = after_long, aes(x = val), colour =  "blue") +
      scale_y_continuous(trans='log2') +
      scale_colour_manual(name = '',
          values =c('orange'='orange','blue'='blue'),
          labels = c('before','after')) +
      xlab(paste(val))


  ggsave(plot = p, filename = paste0(opt$out_dir,"/",filename))


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

#############################################
# Read in mesa files from different sources
#############################################

# Import and combine source metadata files
ls_mesa_meta = apply(df_manifest, 1, importMetaMESA)

# Split into ls of mesa and metadata files for each data set
ls_mesa_inc_count_files <- ls_mesa_allPS_files <- ls_mesa_cluster_files <- c()
ls_mesa_IR_table_files <- ls_mesa_IR_cov_dir <- ls_meta <- ls_sample_names <- c()
for (item in ls_mesa_meta) {
     ls_mesa_inc_count_files <- append(ls_mesa_inc_count_files, item[1])
     ls_mesa_allPS_files <- append(ls_mesa_allPS_files, item[2])
     ls_mesa_cluster_files <- append(ls_mesa_cluster_files, item[3])
     ls_mesa_IR_table_files <- append(ls_mesa_IR_table_files, item[4])
     ls_mesa_IR_cov_dir <- append(ls_mesa_IR_cov_dir, item[5])
     ls_meta <- append(ls_meta, item[6])
     ls_sample_names <- append(ls_sample_names, item[7])
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

# ####################################################
# # Read in all individual mesa IR files and merge
# ####################################################
# Get list of all _intron_coverage.txt files from all batches
ls_mesa_IR_cov_files <- ls_mesa_IR_cov_file_names <- list()
for (dir in ls_mesa_IR_cov_dir){
  files <- file.path(dir, list.files(dir))
  names <- unlist(strsplit(as.character(strsplit(list.files(dir),split= "/")),
                 split = "_intron_coverage.txt"))

  ls_mesa_IR_cov_files <- append(ls_mesa_IR_cov_files,files )
  ls_mesa_IR_cov_file_names <- append(ls_mesa_IR_cov_file_names,names )
}

print(length(ls_mesa_IR_cov_files))
names(ls_mesa_IR_cov_files) <- ls_mesa_IR_cov_file_names

# Files to dfs
ls_mesa_IR_cov_dfs <- foreach(i=ls_mesa_IR_cov_files,
                            n = names(ls_mesa_IR_cov_files), .packages=  'magrittr' ) %dopar% {
    df <- readr::read_delim(i, delim = "\t",   col_names = FALSE) %>%
        dplyr::mutate(X1  = as.character(X1),
                      X2  = as.character(X2),
                      X3  = as.character(X3),
                      X4  = as.character(X4)  ) %>%
        dplyr::select( "X1", "X2","X3","X4", "X5", "X6") %>%
        dplyr::select(-X5,X5) # Move data col to last col

    names(df)[6] <- n # Rename col to sample name
    return(df)
    }

print(length(ls_mesa_IR_cov_dfs))

# Combine all files by columns 
df_merged_IR_cov <- ls_mesa_IR_cov_dfs %>% 
      purrr::reduce(dplyr::inner_join,
                    by = c("X1","X2","X3","X4","X6"),
                    na_matches = "never")  

print(dim(df_merged_IR_cov))

# Write merged table file 
write.table(
  df_merged_IR_cov,
  file.path(opt$out_dir,"merged_mesa_intron_coverage.tsv"), quote=F,sep="\t",
     na="nan", col.names = NA, row.names= TRUE)

# Write back into individual files ; each data col to a different file 
df_bed_cols <- df_merged_IR_cov[, c(1:5)]
ls_cols <- 6:ncol(df_merged_IR_cov)

dir.create(file.path(opt$out_dir,"/mesa_intron_coverage/"))

foreach(i=ls_cols ) %dopar% {
  str_sample_id <- colnames(df_merged_IR_cov)[i]

  # Bind bed columns (same for all) to this col from loop 
  df_IR <- cbind(df_bed_cols, df_merged_IR_cov[,i]) %>%
    dplyr::select(-X6,X6) # Move strand col to last col

  # Add empty cols that may be expected but arent needed when not using mesa table -r
  df_IR$X7 <- df_IR$X8 <- paste0("0,0,0,0,0")

  # Write to file
  write.table(
    x = df_IR,
    file = paste0(opt$out_dir,"/mesa_intron_coverage/",
      str_sample_id,"_intron_coverage.txt"),
    sep="\t",quote=F, col.names = FALSE, row.names = FALSE)

 }

# #########################################
# # Merge mesa inclusion count files
# #########################################
# Note: Do not merge PS files or IR table files because
# they need to be recalculated from merged counts/coverage
mesa_file_names <- list(
              list(ls_mesa_inc_count_files,
                  "manifest_mesa_inclusionCounts_files.tsv","merged_mesa_inclusionCounts.tsv" ))
              # list(ls_mesa_IR_table_files,
              #    "manifest_mesa_IR_files.tsv","merged_mesa_ir_table_intron_retention.tsv" ))

merge <- foreach(i=mesa_file_names ) %dopar% {

        # Write manifest needed for mesa select
        write.table(as.data.frame(unlist(i[1])),
          file=paste0(opt$out_dir,"/",i[2]),
          col.names = FALSE,
          row.names = FALSE,
          quote = FALSE)

        # Run mesa select to merge the mesa file
        cmd <- paste0(
          "mesa select -m ",
                  opt$out_dir,"/",i[2], " -o ",
                  opt$out_dir,"/",i[3], "  --join intersection 2>&1"
                )

        print(cmd)
        system(cmd)
    }

# # Import each merged MESA file and check they have same number of samples
df_merged_inc_counts <- read.table(paste0(opt$out_dir,"/merged_mesa_inclusionCounts.tsv"),
                                    row.names = 1, header=T)
# ######################################
# # Convert Merged count file to PS
# ######################################
# 2>&1 sends standard error standard output
PS_cmd <- paste0(
    "mesa counts_to_ps -i ",
    opt$out_dir,"/merged_mesa_inclusionCounts.tsv -o ",
    opt$out_dir,"/merged_mesa  --recluster 2>&1"
    )
print(PS_cmd)
system(PS_cmd)

# Read in new PS file calculated from merged inclusion counts
df_merged_allPS <- read.table(paste0(opt$out_dir,"/merged_mesa_allPS.tsv"),
                              row.names = 1, header=T)

# #####################################################
# #  Convert Merged non corrected IR cov file to IR PS
# #####################################################
# Run IR table 
cmd <- paste0("mesa ir_table -i ", opt$out_dir, 
        "/merged_mesa_inclusionCounts.tsv -c ",
        opt$out_dir,"/merged_mesa_allClusters.tsv -d ",
        opt$out_dir, "/mesa_intron_coverage -o ",
        opt$out_dir,"/merged_mesa_ir_table")


print(cmd)
system(cmd)
df_merged_IR_PS <- read.table(paste0(opt$out_dir,"/merged_mesa_ir_table_intron_retention.tsv"),
                              row.names = 1, header=T) 
print("df_merged_IR_PS")
print(head(df_merged_IR_PS))
print(dim(df_merged_IR_PS))


if(all.equal(ncol(df_merged_allPS),
            ncol(df_merged_inc_counts),
            ncol(df_merged_IR_PS)) != TRUE)
            stop("Error: Number of columns(samples) in merged inclusion counts, PS, and intron retention are not equal")

# ######################################
# # UMAPs before batch correction
# ######################################
df_to_UMAP(df_merged_inc_counts, "UMAPs_pre_batch_correction/inclusionCounts_PCA_UMAP")
df_to_UMAP(df_merged_allPS, "UMAPs_pre_batch_correction/allPS_PCA_UMAP")
df_to_UMAP(df_merged_IR_PS, "UMAPs_pre_batch_correction/IR_PCA_UMAP")

# # ########################################################
# # # Batch correction of IR coverage + write files
# # #########################################################
dir.create(file.path(opt$out_dir,"/mesa_intron_coverage_bc/"))
              
df_merged_IR_cov_counts <- read.table(paste0(opt$out_dir,
                                  "/merged_mesa_intron_coverage.tsv"),
                                  row.names = 1, header=T)
vals_merged_IR_cov_counts <- df_merged_IR_cov_counts[,-c(1:5)]
# Log2 + 1 transform counts
df_log2trans_IR_cov <- as.data.frame(log2(vals_merged_IR_cov_counts +1))

print(head(df_log2trans_IR_cov))
print(dim(df_log2trans_IR_cov))


# Limma remove batch effect
df_mesa_IR_cov_merge_log2_batch_corr <- limma::removeBatchEffect(
                                  df_log2trans_IR_cov,
                                  batch = df_merged_metadata$data_source,
                                  batch2 = df_merged_metadata$type
                                  )

print(dim(df_mesa_IR_cov_merge_log2_batch_corr))

# Undo log2(x+1) with 2^x - 1
df_mesa_IR_cov_merge_log2_batch_corr = 2^df_mesa_IR_cov_merge_log2_batch_corr -1

# Round all counts below 1 to 0
df_mesa_IR_cov_merge_log2_batch_corr[df_mesa_IR_cov_merge_log2_batch_corr < 1 ] <- 0

df_mesa_IR_cov_merge_log2_batch_corr <- as.data.frame(df_mesa_IR_cov_merge_log2_batch_corr)
print(dim(df_mesa_IR_cov_merge_log2_batch_corr))

head(nrow(df_merged_IR_cov_counts))
head(nrow(df_mesa_IR_cov_merge_log2_batch_corr))

stopifnot(nrow(df_merged_IR_cov_counts)==nrow(df_mesa_IR_cov_merge_log2_batch_corr))


if(nrow(df_merged_IR_cov_counts)!=nrow(df_mesa_IR_cov_merge_log2_batch_corr)) 
  stop("Error: Number of rows different before and after IR batch correction")

# Write Full table 
write.table(
  df_mesa_IR_cov_merge_log2_batch_corr,
  file.path(opt$out_dir,"batch_corr_mesa_intron_coverage.tsv"),
  sep="\t",quote=F, col.names = NA, row.names = TRUE)

# Write back into individual files ; each data col to a different file 
df_bed_cols <- df_merged_IR_cov_counts[, c(1:5)]
ls_cols <- 1:ncol(df_mesa_IR_cov_merge_log2_batch_corr)

print(dim(df_bed_cols))
print(length(ls_cols))
# Loop over IR cov columns
foreach (i = ls_cols , .packages=  'magrittr') %dopar% {      
  str_sample_id <- colnames(df_mesa_IR_cov_merge_log2_batch_corr)[i]

  # Bind bed columns (same for all) to this col from loop 
  df_IR <- cbind(df_bed_cols, df_mesa_IR_cov_merge_log2_batch_corr[,i]) %>%
    dplyr::select(-X6,X6) # Move strand col to last col

  # Add empty rows that maybe expected but arent needed
  df_IR$X7 <- df_IR$X8 <- paste0("0,0,0,0,0")

  # Write to file
  write.table(
    x = df_IR,
    file = paste0(opt$out_dir,"/mesa_intron_coverage_bc/",
      str_sample_id,"_intron_coverage.txt"),
    sep="\t",quote=F, col.names = FALSE, row.names = FALSE)
}


# ########################################################
# # Batch correction of inclusion counts
# #########################################################

# Log2 + 1 transform counts
df_log2trans_inc_counts <- as.data.frame(log2(df_merged_inc_counts +1))

# Limma remove batch effect
df_mesa_inc_count_merge_log2_batch_corr <- limma::removeBatchEffect(
                                  df_log2trans_inc_counts,
                                  batch = df_merged_metadata$data_source,
                                  batch2 = df_merged_metadata$type
                                  )

# Undo log2(x+1) with 2^x - 1
df_merged_inc_counts_batch_corr = 2^df_mesa_inc_count_merge_log2_batch_corr -1
rownames(df_merged_inc_counts_batch_corr) <- rownames(df_merged_inc_counts)

# Round all counts below 1 to 0
df_merged_inc_counts_batch_corr[df_merged_inc_counts_batch_corr < 1 ] <- 0

print(head(df_merged_inc_counts_batch_corr))
print(dim(df_merged_inc_counts_batch_corr))


df_merged_inc_counts_batch_corr <- as.data.frame(df_merged_inc_counts_batch_corr) %>% 
  tibble::rownames_to_column("cluster")


write.table(
  df_merged_inc_counts_batch_corr,
  file.path(opt$out_dir,"batch_corr_mesa_inclusionCounts.tsv"),
  sep="\t",quote=F, row.names = FALSE)

df_merged_inc_counts_batch_corr <- as.data.frame(df_merged_inc_counts_batch_corr) %>% 
  tibble::column_to_rownames("cluster")

# ########################################################
# # Convert batch corrected counts to PS
# #########################################################
# # Merge all PS files with MESA select_samples command
# # 2>&1 sends standard error standard output
PS_cmd <- paste0(
    "mesa counts_to_ps -i ",
    opt$out_dir,"/batch_corr_mesa_inclusionCounts.tsv -o ",
    opt$out_dir,"/batch_corr_mesa  --recluster 2>&1"
    )
print(PS_cmd)
system(PS_cmd)
########################################################
# Convert batch corrected IR counts to IR PS
#########################################################
# Run IR table to convert from coverage to ???
# Run mesa IR table to merge from counts to IR 

cmd <- paste0("mesa ir_table -i ", opt$out_dir, 
        "/batch_corr_mesa_inclusionCounts.tsv -c ",
        opt$out_dir,"/merged_mesa_allClusters.tsv -d ",
        opt$out_dir, "/mesa_intron_coverage_bc -o ",
        opt$out_dir,"/batch_corr_mesa_ir_table  2>&1")

print(cmd)
system(cmd)


df_merged_IR_batch_corr <- read.table(paste0(opt$out_dir,
                              "/batch_corr_mesa_ir_table_intron_retention.tsv"),
                              row.names = 1, header=T) 

print(head(df_merged_IR_batch_corr))
print(dim(df_merged_IR_batch_corr))
########################################################
# Open PS after batch correction / write LM22 subset 
#########################################################
df_merged_allPS_batch_corr <- read.table(paste0(opt$out_dir,"/batch_corr_mesa_allPS.tsv"),
                              row.names = 1, header=T)

# Drop non LM22 samples from mesa PS
df_merged_allPS_batch_corr_lm22 <- df_merged_allPS_batch_corr %>%
  dplyr::select(ls_smpls_lm22)
rownames(df_merged_allPS_batch_corr_lm22) <- rownames(df_merged_allPS_batch_corr)

write.table(
  df_merged_allPS_batch_corr_lm22,
  file.path(opt$out_dir,"batch_corr_mesa_allPS_LM22.tsv"), quote=F,sep="\t", na="nan",
  col.names = NA, row.names= TRUE)

######################################
# UMAPs after batch correction
######################################
df_to_UMAP(df_merged_inc_counts_batch_corr, "UMAPs_post_batch_correction/inclusionCounts_PCA_UMAP")
df_to_UMAP(df_merged_allPS_batch_corr, "UMAPs_post_batch_correction/allPS_PCA_UMAP")
df_to_UMAP(df_merged_IR_batch_corr, "UMAPs_post_batch_correction/IR_PCA_UMAP")

#################################################
# Compare distributions before and after BC
#################################################
plot_before_after(df_merged_inc_counts, as.data.frame(df_merged_inc_counts_batch_corr), "hist_inclusionCounts_before_after_bc.png", "Inclusion count")
plot_before_after(df_merged_allPS, as.data.frame(df_merged_allPS_batch_corr), "hist_PS_before_after_bc.png","PS")
plot_before_after(df_merged_IR_PS, as.data.frame(df_merged_IR_batch_corr), "hist_IR_before_after_bc.png","Intron coverage")