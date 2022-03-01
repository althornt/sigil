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
  res_IR_path <- file.path(row[2], "mesa_out", "mesa_ir_table_intron_retention.tsv")


  return(list(
      "ls_mesa_inc_count_files"=res_inc_count_path,
      "ls_mesa_allPS_files"=res_allPS_path,
      "ls_mesa_cluster_files"=res_cluster_path,
      "ls_mesa_IR_files" = res_IR_path,
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

plot_before_after <- function(df_before, df_after, filename, val){

  before_long <- df_before %>%
  # dplyr::filter(rownames(df_before) %in% c("21:33432871-33436827:+","21:33432871-33436827:-")) %>%
    tibble::rownames_to_column("event") %>%
    tidyr::gather(., key="sample", value = val, -c(event)) %>%
    as.data.frame()

  after_long <- df_after %>%
    # dplyr::filter(rownames(df_after) %in% c("21:33432871-33436827:+","21:33432871-33436827:-")) %>%
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


  ggsave(plot = p, filename = paste0(opt$out_dir,"/",filename,"_IFNGR2.png"))


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
ls_mesa_IR_files <- ls_meta <- ls_sample_names <- c()
for (item in ls_mesa_meta) {
     ls_mesa_inc_count_files <- append(ls_mesa_inc_count_files, item[1])
     ls_mesa_allPS_files <- append(ls_mesa_allPS_files, item[2])
     ls_mesa_cluster_files <- append(ls_mesa_cluster_files, item[3])
     ls_mesa_IR_files <- append(ls_mesa_IR_files, item[4])
     ls_meta <- append(ls_meta, item[5])
     ls_sample_names <- append(ls_sample_names, item[6])
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

#########################################
# Merge mesa files
#########################################

mesa_file_names <- list(
              list(ls_mesa_inc_count_files,
                  "manifest_mesa_inclusionCounts_files.tsv","merged_mesa_inclusionCounts.tsv" ),
              list(ls_mesa_allPS_files,
                 "manifest_mesa_allPS_files.tsv","merged_mesa_allPS.tsv" ),
              list(ls_mesa_IR_files,
                 "manifest_mesa_IR_files.tsv","merged_mesa_ir_table_intron_retention.tsv" ))

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
        # system(cmd)
    }

# Import each merged MESA file and check they have same number of samples
df_merged_inc_counts <- read.table(paste0(opt$out_dir,"/merged_mesa_inclusionCounts.tsv"),
                                    row.names = 1, header=T)
df_merged_allPS <- read.table(paste0(opt$out_dir,"/merged_mesa_allPS.tsv"),
                              row.names = 1, header=T)
df_merged_IR <- read.table(paste0(opt$out_dir,"/merged_mesa_ir_table_intron_retention.tsv"),
                                    row.names = 1, header=T)
if(all.equal(ncol(df_merged_allPS),
            ncol(df_merged_inc_counts),
            ncol(df_merged_IR)) != TRUE)
  stop("Error: Number of columns across PS, inclusion counts, and intron retention files are not equal")

######################################
# UMAPs before batch correction
######################################
df_to_UMAP(df_merged_inc_counts, "UMAPs_pre_batch_correction/inclusionCounts_PCA_UMAP")
df_to_UMAP(df_merged_allPS, "UMAPs_pre_batch_correction/allPS_PCA_UMAP")
df_to_UMAP(df_merged_IR, "UMAPs_pre_batch_correction/IR_PCA_UMAP")

########################################################
# Batch correction of inclusion counts
#########################################################
# Log2 + 1 transform counts
df_log2trans_inc_counts <- as.data.frame(log2(df_merged_inc_counts +1))
plot_PS_hist(df_merged_inc_counts,  paste0(opt$out_dir,"/hist_merged_inc_counts.png"))
plot_PS_hist(df_log2trans_inc_counts,  paste0(opt$out_dir,"/hist_merged_inc_counts_log.png"))

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

plot_PS_hist(as.data.frame(df_merged_inc_counts_batch_corr),
              paste0(opt$out_dir,"/hist_merged_inc_counts_bc.png"))

write.table(
  df_merged_inc_counts_batch_corr,
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
    opt$out_dir,"/batch_corr_mesa  --recluster 2>&1"
    )
print(PS_cmd)
system(PS_cmd)


########################################################
# UMAPs of PS after batch correction
#########################################################
df_merged_allPS_batch_corr <- read.table(paste0(opt$out_dir,"/batch_corr_mesa_allPS.tsv"),
                              row.names = 1, header=T)

plot_PS_hist(df_merged_allPS_batch_corr,  paste0(opt$out_dir,"/hist_allPS_bc.png"))

# Drop non LM22 samples from mesa PS
df_merged_allPS_batch_corr_lm22 <- df_merged_allPS_batch_corr %>%
  dplyr::select(ls_smpls_lm22)
rownames(df_merged_allPS_batch_corr_lm22) <- rownames(df_merged_allPS_batch_corr)

write.table(
  df_merged_allPS_batch_corr_lm22,
  file.path(opt$out_dir,"batch_corr_mesa_allPS_LM22.tsv"), quote=F,sep="\t", na="nan",
  col.names = NA, row.names= TRUE)

######################################
# Intron Retention batch correction
######################################
# Log2 + 1 transform IR
df_merged_IR_log2 <- as.data.frame(log2(df_merged_IR +1))

# Limma remove batch effect
df_merged_IR_log2_batch_corr <- limma::removeBatchEffect(
                                  df_merged_IR_log2,
                                  batch = df_merged_metadata$data_source,
                                  batch2 = df_merged_metadata$type
                                  )

# Undo log2(x+1) with 2^x - 1
df_merged_IR_batch_corr = 2^df_merged_IR_log2_batch_corr -1
rownames(df_merged_IR_batch_corr) <- rownames(df_merged_IR)

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
plot_before_after(df_merged_IR, as.data.frame(df_merged_IR_batch_corr), "hist_IR_before_after_bc.png","Intron coverage")

#################################
# IFNGR2 example
################################

# inc_count_long <- df_merged_inc_counts %>%
#   dplyr::filter(rownames(df_merged_inc_counts) %in% c("21:33432871-33436827:+","21:33432871-33436827:-")) %>%
#   tibble::rownames_to_column("event") %>%
#   tidyr::gather(., key="sample", value = "Inc_counts", -c(event)) %>%
#   as.data.frame()
#
# bc_inc_count_long <- data.frame(df_merged_inc_counts_batch_corr) %>%
#   dplyr::filter(rownames(df_merged_inc_counts_batch_corr) %in% c("21:33432871-33436827:+","21:33432871-33436827:-")) %>%
#   tibble::rownames_to_column("event") %>%
#   tidyr::gather(., key="sample", value = "Inc_counts", -c(event)) %>%
#   as.data.frame()
#
#
# p <- ggplot() +
#     geom_freqpoly(data = inc_count_long, aes(x = Inc_counts), colour =  "orange") +
#     geom_freqpoly(data = bc_inc_count_long, aes(x = Inc_counts), colour =  "blue") +
#     scale_y_continuous(trans='log2') +
#     scale_colour_manual(name = '',
#         values =c('orange'='orange','blue'='blue'),
#         labels = c('before','after'))
#
# ggsave(plot = p, filename = paste0(opt$out_dir,"/hist_PS_before_after_bc_ifngr2_count.png"))
