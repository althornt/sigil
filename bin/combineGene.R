#!/usr/bin/env Rscript

library(optparse)
library(magrittr)
library(dplyr)
library(uwot)
library(ggplot2)
library(RColorBrewer)
library(ensembldb)
library(DESeq2)
library(pheatmap)

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
    help = "path to tsv manifest file with 3 columns: data_source, res_path, metadata_path")

  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directory
if (!dir.exists(file.path(opt$out_dir,"/UMAPs/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs/"), recursive = TRUE, showWarnings = TRUE)}

# Open files
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(head(df_manifest))

# Check all combined kallisto results files exist
ls_kallisto_paths <- file.path(
  df_manifest$res_path,
  "post_kallisto_out",
   "combined_kallisto_log2tpm.csv")
if(all(file.exists(ls_kallisto_paths)) != TRUE)
  stop("Error: missing combined_kallisto_log2tpm.csv files.")

importKallistoMeta <- function(row){
  source_name <- row[1]

  # Read kallisto file
  res_path <- file.path(row[2], "post_kallisto_out","combined_kallisto_log2tpm.csv")
  df_res <- read.csv(res_path) %>%
    dplyr::rename("gene" = "X")

  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general, LM22)

  df_metadata$data_source <- source_name
  print(df_metadata)

  return(list("res"=df_res, "metadata"=df_metadata))
}

# Import and combine kallisto files
ls_kallisto_meta = apply(df_manifest, 1, importKallistoMeta)

# Split into kallisto and metadata files for each data set
ls_kallisto <- ls_mesa <- c()
for (item in ls_kallisto_meta) {
     ls_kallisto <- append(ls_kallisto, item[1])
     ls_mesa <- append(ls_mesa, item[2])
   }

# Merge metadata by rows
df_merged_metadata <- do.call("rbind", ls_mesa)
rownames(df_merged_metadata) <- c()

# Merge kallisto by columns
df_merged_kallisto_res <- Reduce(function(...) merge(..., all = TRUE, by="gene"), ls_kallisto)
rownames(df_merged_kallisto_res) <- df_merged_kallisto_res$gene

# Convert to numeric
df_merged_kallisto_res <- dplyr::mutate_all(
    df_merged_kallisto_res,
    function(x) as.numeric(as.character(x)))
write.csv(df_merged_kallisto_res,
          file.path(opt$out_dir,"combined_kallisto_log2tpm.csv"),
          row.names = TRUE)

# Make dfs for only samples with LM22 labels
ls_smpls_lm22 <- df_merged_metadata %>%
  dplyr::filter(LM22 != "") %>%
  dplyr::pull(Run)
df_merged_metadata_lm22 <- df_merged_metadata %>%
  dplyr::filter(LM22 != "")
df_merged_kallisto_res_lm22 <- df_merged_kallisto_res[ls_smpls_lm22]
write.csv(df_merged_kallisto_res_lm22,
          file.path(opt$out_dir,"combined_kallisto_log2tpm_lm22.csv"),
          row.names = TRUE)

# UMAP function
make_umap <- function(num_neighbor,meta_col) {

  set.seed(123)

  # Make color palette
  n <- length(unique(df_merged_metadata_lm22[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, df_merged_metadata_lm22)

  # Plot UMAP
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP Kallisto: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste("UMAPs/kallisto_PCA_UMAP",meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)

}

# Drop genes with low variance.
getVar <- apply(df_merged_kallisto_res_lm22[, -1], 1, var)
param <- median(getVar)
log2trans_dat_filt <- df_merged_kallisto_res_lm22[getVar > param & !is.na(getVar), ]

# Transpose and format
log2trans_dat_filt_t <- as.data.frame(t(log2trans_dat_filt)) %>%
  dplyr::filter(!row.names(.) %in% c("gene"))
rownames(log2trans_dat_filt_t) <- colnames(log2trans_dat_filt)[-1] # skip gene

# PCA.
prcomp.out = as.data.frame(prcomp(log2trans_dat_filt_t, scale = F)$x)
prcomp.out$Run = rownames(log2trans_dat_filt_t)
prcomp.out.merge = merge(prcomp.out, y = df_merged_metadata_lm22)

# Making variations of UMAPs with different numbers of neighbors
lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22")
lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general")
lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source")
