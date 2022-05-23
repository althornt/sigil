#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(uwot)
library(ggplot2)
library(RColorBrewer)
# library(foreach)
# library(doParallel)
# cl <- makeCluster(detectCores() - 1, outfile = "")
# registerDoParallel(cl)

importMetaMESA <- function(row){
  # Read metadata and add column for run
  df_metadata <- read.csv(file.path(row[3])) %>%
    dplyr::select(Run, sigil_general, main_label, group_label)

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
      "ls_samples_run"=df_metadata$Run,
      "source"=unique(df_metadata$data_source)))
}

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
    help = "path to metadata file"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Open manifest
df_manifest <- read.csv(file = opt$manifest, sep = "\t", header=TRUE)
print(df_manifest)

dir.create(file.path(opt$out_dir), recursive = TRUE, showWarnings = FALSE)


# Import and combine source metadata files
ls_mesa_meta = apply(df_manifest, 1, importMetaMESA)

# Split into ls of mesa and metadata files for each data set
ls_mesa_inc_count_files <- ls_mesa_allPS_files <- ls_mesa_cluster_files <- c()
ls_mesa_IR_table_files <- ls_mesa_IR_cov_dir <- ls_meta <- ls_sample_names <- c()
ls_source_names <- c()
for (item in ls_mesa_meta) {
     ls_mesa_inc_count_files <- append(ls_mesa_inc_count_files, item[1])
     ls_mesa_allPS_files <- append(ls_mesa_allPS_files, item[2])
     ls_mesa_cluster_files <- append(ls_mesa_cluster_files, item[3])
     ls_mesa_IR_table_files <- append(ls_mesa_IR_table_files, item[4])
     ls_mesa_IR_cov_dir <- append(ls_mesa_IR_cov_dir, item[5])
     ls_meta <- append(ls_meta, item[6])
     ls_sample_names <- append(ls_sample_names, item[7])
     ls_source_names <- append(ls_source_names, item[8]$source)

   }

# Combine metadata from each data source by rows
df_merged_metadata <- do.call("rbind", ls_meta)
rownames(df_merged_metadata) <- c()

# print(ls_source_names)
# print("--")

# print(ls_source_names[3])

# names(ls_mesa_allPS_files) <- ls_source_names
# print(ls_mesa_allPS_files)

# print(names(ls_mesa_allPS_files))

# quit()

# columns= c("PS_nrow","PS_ncol") 
# df_summary = data.frame(matrix( ncol = length(columns))) 
# assign column names
# colnames(df_summary) = columns
# rownames(df_summary) = ls_source_names

df_summary = data.frame() 
print(df_summary)

for (i in seq_along(ls_mesa_allPS_files)){
    source_name = ls_source_names[i]
    file_name = ls_mesa_allPS_files[i]
    print(source_name)

    df_allPS <- read.table(as.character(file_name), row.names = 1, header=T)
    print(dim(df_allPS))

    df_summary[source_name,"PS_nrow" ] <- nrow(df_allPS)
    df_summary[source_name,"PS_ncol" ] <- ncol(df_allPS)

    #count nan per sample and per junction
    column.nan.counts <- colSums(is.na(df_allPS[,2:ncol(df_allPS)]))
    row.nan.counts <- rowSums(is.na(df_allPS))
    high_nan_samples = names(column.nan.counts[column.nan.counts > (dim(df_allPS)[1]/2)])

    print(min(column.nan.counts))
    print(max(column.nan.counts))
    print(median(column.nan.counts))

    df_summary[source_name,"PS_min_nan_col_count" ] <- min(column.nan.counts)
    df_summary[source_name,"PS_max_nan_col_count" ] <- max(column.nan.counts)
    df_summary[source_name,"PS_med_nan_col_count" ] <- median(column.nan.counts)

    # Percent of nan in all df 
    df_summary[source_name,"PS_perc_nan" ]<-sum(is.na(df_allPS))/prod(dim(df_allPS))
    
    df_summary[source_name,"PS_min_nan_col_perc" ] <- min(colMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_max_nan_col_perc" ] <- max(colMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_med_nan_col_perc" ] <- median(colMeans(is.na(df_allPS)))

    #plot dist of nan per column
    ggplot(NULL, aes(x=colMeans(is.na(df_allPS)))) + geom_histogram() + 
            scale_x_continuous(limits = c( 0,1 )) + xlab("Percent NAN in sample") +
            ggtitle(source_name) + geom_vline(xintercept = 0.5)

    ggsave(file.path(opt$out_dir,
            paste(source_name,"sample_nans", "png", sep = '.')),
            device = "png",
            width = 12,
            dpi = 300)

    # Percent nan in each junction
    df_summary[source_name,"PS_min_nan_row_perc" ] <- min(rowMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_max_nan_row_perc" ] <- max(rowMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_med_nan_row_perc" ] <- median(rowMeans(is.na(df_allPS)))

    #plot dist of nan per junction
    ggplot(NULL, aes(x=rowMeans(is.na(df_allPS)))) + geom_histogram() + 
            scale_x_continuous(limits = c( 0,1 )) + xlab("Percent NAN in junction") +
            ggtitle(source_name) + geom_vline(xintercept = 0.5)
    ggsave(file.path(opt$out_dir,
            paste(source_name,"junction_nans", "png", sep = '.')),
            device = "png",
            width = 12,
            dpi = 300)


}

print(df_summary)