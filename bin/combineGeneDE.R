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
library(tidyverse)
library(limma)
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

##################
# Functions
##################

importMetaKallisto <- function(row){
  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general, LM22, LM6)

  # Add metadata to column
  df_metadata$data_source <- row[1] # add name of data source
  df_metadata$type <- row[4] # add rna-seq type (paired vs single)

  # Get paths to kallisto files
  res_path <- file.path(row[2], "kallisto_out",df_metadata$Run, "abundance.h5")

  return(list(
      "ls_kallsto_paths"=res_path,
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
    labs(title= paste("UMAP Kallisto: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste(out_path,meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)

}

# Run Deseq2 for one cell type vs all other samples
runDE_1_vs_all <- function(meta_col_to_use, cell_type_val) {
  # Make sample table to compare given cell type vs all other cell types
  sampleTable <- df_merged_metadata_lm22 %>%
    dplyr::select(Run, meta_col_to_use, data_source, type) %>%
    dplyr::mutate(condition = ifelse(get(meta_col_to_use) == cell_type_val, "main", "other"))

  path_to_deseq2_folder <- paste0(opt$out_dir,"/compare_",meta_col_to_use,"/deseq2_outputs/")
  runDEseq2(sampleTable,meta_col_to_use, cell_type_val,path_to_deseq2_folder)
}

runDEseq2 <- function(sampleTable,meta_col_to_use, cell_type_val,path_to_deseq2_folder){

  if (!dir.exists(file.path(path_to_deseq2_folder))){
    dir.create(file.path(path_to_deseq2_folder),recursive = TRUE, showWarnings = TRUE)
    }

  #Write output to file for checking (not read by script)
  write.csv(sampleTable,file.path(paste(path_to_deseq2_folder,cell_type_val,"sampletable.csv")))
  print(sampleTable)

  ls_kallisto_paths_lm22_ <- ls_kallisto_paths_lm22[rownames(sampleTable)]

  # If all data_source and seq type the same 
  # If all data_source the same 
  # If all seq type the same

  print(unique(sampleTable$data_source))
  print(unique(sampleTable$type))


  txi.kallisto <- tximport(ls_kallisto_paths_lm22_, type = "kallisto",tx2gene = tx2gene,
                    txOut = FALSE,ignoreTxVersion=TRUE,
                    countsFromAbundance = "scaledTPM")

  # Run DEseq and base the model sequencing type and  data source if possible
  if ((length(unique(sampleTable$data_source))==1) & (length(unique(sampleTable$type))==1)){
      print(cell_type_val)
      print("only condition")
      dds <- DESeqDataSetFromTximport(
              txi = txi.kallisto, 
              colData = sampleTable,
              design =  ~condition)
  }  else if ((length(unique(sampleTable$data_source))==1) & (length(unique(sampleTable$type))>1)){
      print(cell_type_val)
      print("only condition, type")
      dds <-  tryCatch({
        dDESeqDataSetFromTximport(
                    txi = txi.kallisto, 
                    colData = sampleTable,
                    design =  ~type + condition)
      },error = function(err){
        DESeqDataSetFromTximport(
                    txi = txi.kallisto, 
                    colData = sampleTable,
                    design =  ~ condition)
      })
  }  else if ((length(unique(sampleTable$data_source))>1) & (length(unique(sampleTable$type))==1)){
      print(cell_type_val)
      print("only condition, data_source")
      dds <-  tryCatch({
        dDESeqDataSetFromTximport(
                    txi = txi.kallisto, 
                    colData = sampleTable,
                    design =  ~ data_source + condition)
      },error = function(err){
        DESeqDataSetFromTximport(
                    txi = txi.kallisto, 
                    colData = sampleTable,
                    design =  ~ condition)
      })
  }  else if ((length(unique(sampleTable$data_source))>1) & (length(unique(sampleTable$type))>1)){
      print(cell_type_val)
      print("condition, type, and data_source")
      dds <-  tryCatch({
        dDESeqDataSetFromTximport(
                    txi = txi.kallisto, 
                    colData = sampleTable,
                    design =  ~ data_source + type + condition)
      },error = function(err){
        DESeqDataSetFromTximport(
                    txi = txi.kallisto, 
                    colData = sampleTable,
                    design =  ~ condition)
      })
  }
  
  # run DEseq2 on main cell types vs other
  colData(dds)$condition<-factor(
    colData(dds)$condition,
    levels=c("other", "main"))
  dds_<-DESeq(dds)
  res<-results(dds_)
  res<-res[order(res$padj),]

  #Write output to file
  write.csv(as.data.frame(res),file.path(paste(path_to_deseq2_folder,cell_type_val,".csv")))
}

runDE_1_vs_all_within_type <- function(ls_cell_types, label) {
  ls_cell_types <- unlist(ls_cell_types)  
  print(ls_cell_types)
  for (cell_type_val in ls_cell_types){
    # Make sample table to compare given cell type vs all other cell types
    sampleTable <- df_merged_metadata_lm22 %>%
      dplyr::filter(LM22 %in% ls_cell_types) %>%
      dplyr::select(Run, LM22, data_source, type) %>%
      dplyr::mutate(condition = ifelse(LM22 == cell_type_val, "main", "other")) %>%
      tibble::column_to_rownames("Run")
    path_to_deseq2_folder <- paste0(opt$out_dir,"/compare_within_type/deseq2_outputs/")
    runDEseq2(sampleTable,"LM22", cell_type_val,path_to_deseq2_folder)
  }  
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


#####################
# Processing inputs
#####################

# Open files
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(head(df_manifest))

# Import and combine source metadata files
ls_kallisto_meta = apply(df_manifest, 1, importMetaKallisto)

# Split into kallisto and metadata files for each data set
ls_kallisto <- ls_meta <- ls_sample_names <- c()
for (item in ls_kallisto_meta) {
     ls_kallisto <- append(ls_kallisto, item[1])
     ls_meta <- append(ls_meta, item[2])
     ls_sample_names <- append(ls_sample_names, item[3])
   }

# Merge metadata from each data source by rows
df_merged_metadata <- do.call("rbind", ls_meta)
rownames(df_merged_metadata) <- c()

# List of samples with LM22 labels
ls_smpls_lm22 <- df_merged_metadata %>%
   dplyr::filter(LM22 != "") %>%
   dplyr::pull(Run)

# Remove samples without LM22 labels
df_merged_metadata_lm22 <- df_merged_metadata %>%
   dplyr::filter(LM22 != "")
write.csv(df_merged_metadata_lm22,
            file.path(file.path(opt$out_dir,"metadata.csv")),
            row.names = FALSE)
# Make list of kallisto paths from each data source
ls_kallisto_paths <- unlist(ls_kallisto,use.names = FALSE,recursive = FALSE)
names(ls_kallisto_paths) <- unlist(ls_sample_names)

# Remove kallisto paths without LM22 labels
ls_kallisto_paths_lm22 <- ls_kallisto_paths[ls_smpls_lm22]

# Import counts with tximport
if(all(file.exists(ls_kallisto_paths_lm22)) != TRUE)
  stop("Error: missing kallisto h5 files.")

# Transcripts to gene, used in tximport
edb <- EnsDb.Hsapiens.v86
tx2gene = transcripts(edb, columns=c("tx_id", "gene_name"),return.type="DataFrame")

txi.kallisto <- tximport(ls_kallisto_paths_lm22, type = "kallisto",tx2gene = tx2gene,
                         txOut = FALSE,ignoreTxVersion=TRUE,
                         countsFromAbundance = "scaledTPM")
dat = txi.kallisto$counts

# Write combined TPM file and log2 TPM
write.csv(dat,
         file.path(opt$out_dir,"combined_kallisto_tpm.csv"),
         row.names = TRUE)
log2trans_dat <- as.data.frame(log2(dat +1))
write.csv(log2trans_dat,
         file.path(file.path(opt$out_dir,"combined_kallisto_log2tpm.csv")),
         row.names = TRUE)

metadata_summary <- df_merged_metadata_lm22 %>%
  dplyr::count(LM22, data_source)
write.csv(metadata_summary,
           file.path(file.path(opt$out_dir,"metadata_summary.csv")),
           row.names = FALSE)

##################
# DESEQ2 LM22
# ##################
# # Run deseq2 on each LM22 cell type vs all others
# if("LM22" %in% colnames(df_merged_metadata_lm22)){
#   ls_lm22_cell_types <-  unique(df_merged_metadata_lm22[["LM22"]])[1]
#   print(ls_lm22_cell_types)
#   foreach(i=ls_lm22_cell_types, .packages=  c('magrittr', 'DESeq2')) %dopar% {
#     runDE_1_vs_all(
#         cell_type_val = i,
#         meta_col_to_use="LM22"
#         )
#       }
# }

##################
# DESEQ2 LM6
##################
# Run deseq2 on each LM6 cell type vs all others
# if("LM6" %in% colnames(df_merged_metadata_lm22)){
#   ls_lm6_cell_types <-  unique(df_merged_metadata_lm22[["LM6"]])[1]
#   print(ls_lm6_cell_types)
#   print("-----")
#   foreach(i=ls_lm6_cell_types, .packages=  c('magrittr', 'DESeq2')) %dopar% {
#     runDE_1_vs_all(
#         cell_type_val = i,
#         meta_col_to_use="LM6"
#         )
#     }
# }

##################
# Within Cell Type 
##################
# Dont yet include types without samples or else nan filtering will break

# T_cell_types <- list(
#   "T cells CD8",
#   "T cells CD4 naive",
#   "T cells CD4 memory resting",
#   "T cells CD4 memory  activated",
#   "T cells follicular helper",
#   "T cells regulatory (Tregs)",
#   "T cells gamma delta")
T_cell_types <- list(
  "T cells CD8",
  "T cells CD4 naive",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta")

# mon_mac_cell_types <- list(
#   "Monocytes",
#   "Macrophages M0",
#   "Macrophages M1",
#   "Macrophages M2")
mon_mac_cell_types <- list(
    "Monocytes",
    "Macrophages M0",
    "Macrophages M1")

B_cell_types <- list(
  "B cells naive",
  "B cells memory")

dendritic_cell_types <- list(
  "Dendritic cells resting",
  "Dendritic cells activated")

# mast_cell_types <- list(
#   "Mast cells resting",
#   "Mast cells activated")

# NK_cell_types <- list(
#   "NK cells resting",
#   "NK cells activated")

# ls_within_cell_types <- list(
#   "T_cell_types" = T_cell_types,
#   "mon_mac_cell_types" = mon_mac_cell_types,
#   "Bcell" = B_cell_types,
#   "Dendritic" = dendritic_cell_types,
#   "Mast" = mast_cell_types,
#   "NK" = NK_cell_types)

ls_within_cell_types <- list(
  list("Tcell", T_cell_types),
  list("Mon_Mac", mon_mac_cell_types),
  list("Bcell", B_cell_types),
  list("Dendritic" ,dendritic_cell_types)
)

# Run Deseq2 for one cell type vs all other samples

foreach(i=ls_within_cell_types, .packages=  c('magrittr', 'DESeq2','tximport')) %dopar% {
  print("foreach")
  runDE_1_vs_all_within_type(
      ls_cell_types = i[2], 
      label= i[1]
      )
  }

####################################
# UMAPs before batch correction
####################################

# # Drop genes with low variance.
# getVar <- apply(log2trans_dat[, -1], 1, var)
# param <- median(getVar)
# log2trans_dat_filt <- log2trans_dat[getVar > param & !is.na(getVar), ]

# # Transpose and format
# log2trans_dat_filt_t <- as.data.frame(t(log2trans_dat_filt)) %>%
#   dplyr::filter(!row.names(.) %in% c("gene"))
# rownames(log2trans_dat_filt_t) <- colnames(log2trans_dat_filt)

# # PCA.
# prcomp.out = as.data.frame(prcomp(log2trans_dat_filt_t, scale = F)$x)
# prcomp.out$Run = rownames(log2trans_dat_filt_t)
# prcomp.out.merge = merge(prcomp.out, y = log2trans_dat)

# # Making variations of UMAPs with different numbers of neighbors
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/kallisto_PCA_UMAP")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/kallisto_PCA_UMAP")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source",
#   df_PCA = prcomp.out, out_path = "UMAPs_pre_batch_correction/kallisto_PCA_UMAP")

# ####################################
# # UMAPs after batch correction
# ####################################

# # Batch correction
# log2trans_dat_batch_corr <- limma::removeBatchEffect(
#                                   log2trans_dat,
#                                   batch = df_merged_metadata_lm22$data_source,
#                                   batch2 = df_merged_metadata_lm22$type
#                                   )

#   write.csv(log2trans_dat_batch_corr,
#            file.path(file.path(opt$out_dir,"combined_kallisto_log2tpm_batch_corrected.csv")),
#            row.names = TRUE)

# # Drop genes with low variance.
# getVar_batch_corr <- apply(log2trans_dat_batch_corr[, -1], 1, var)
# param <- median(getVar_batch_corr)
# log2trans_dat_batch_corr_filt <- log2trans_dat_batch_corr[getVar_batch_corr > param & !is.na(getVar_batch_corr), ]

# # Transpose and format
# log2trans_dat_batch_corr_filt_t <- as.data.frame(t(log2trans_dat_batch_corr_filt)) %>%
#   dplyr::filter(!row.names(.) %in% c("gene"))
# rownames(log2trans_dat_batch_corr_filt_t) <- colnames(log2trans_dat_batch_corr_filt)

# # PCA.
# prcomp.out.batch.corr = as.data.frame(prcomp(log2trans_dat_batch_corr_filt_t, scale = F)$x)
# prcomp.out.batch.corr$Run = rownames(log2trans_dat_batch_corr_filt_t)
# # prcomp.out.merge = merge(prcomp.out, y = log2trans_dat)

# # Making variations of UMAPs with different numbers of neighbors
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22",
#   df_PCA = prcomp.out.batch.corr, out_path = "UMAPs_post_batch_correction/kallisto_PCA_UMAP")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general",
#   df_PCA = prcomp.out.batch.corr, out_path = "UMAPs_post_batch_correction/kallisto_PCA_UMAP")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source",
#   df_PCA = prcomp.out.batch.corr, out_path = "UMAPs_post_batch_correction/kallisto_PCA_UMAP")