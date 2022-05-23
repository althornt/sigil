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
mesatable_qc <-function(ls_mesa_files, out_prefix){
  
  df_summary = data.frame() 
  for (i in seq_along(ls_mesa_files)){


    source_name = ls_source_names[i]
    print(source_name)

    file_name = as.character(ls_mesa_files[i])
    print(file_name)
    if (!file.exists(file_name)){
      next
    }

    df_allPS <- read.table(as.character(file_name), row.names = 1, header=T)
    print(dim(df_allPS))

    # Add dimensions of PS to summary table 
    df_summary[source_name,"PS_nrow" ] <- nrow(df_allPS)
    df_summary[source_name,"PS_ncol" ] <- ncol(df_allPS)

    # # Count nan per sample and per junction
    # high_nan_samples = names(column.nan.counts[column.nan.counts > (dim(df_allPS)[1]/2)])

    # # df_summary[source_name,"PS_min_nan_col_count" ] <- min(column.nan.counts)
    # # df_summary[source_name,"PS_max_nan_col_count" ] <- max(column.nan.counts)
    # # df_summary[source_name,"PS_med_nan_col_count" ] <- median(column.nan.counts)

    # Percent of nan in all df 
    df_summary[source_name,"PS_perc_nan" ]<-sum(is.na(df_allPS))/prod(dim(df_allPS))
    
    # Percent nan in each sample 
    df_summary[source_name,"PS_min_nan_col_perc" ] <- min(colMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_max_nan_col_perc" ] <- max(colMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_med_nan_col_perc" ] <- median(colMeans(is.na(df_allPS)))

    # Plot hist of nan per sample in given data souce 
    ggplot(NULL, aes(x=colMeans(is.na(df_allPS)))) + geom_histogram() + 
            scale_x_continuous(limits = c( 0,1 )) + xlab("Percent NAN in sample") +
            ggtitle(source_name) + geom_vline(xintercept = 0.5)

    ggsave(file.path(opt$out_dir,
            paste(out_prefix, source_name,"sample_nans", "png", sep = '.')),
            device = "png",
            width = 12,
            dpi = 300)

    # Percent nan in each junction
    df_summary[source_name,"PS_min_nan_row_perc" ] <- min(rowMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_max_nan_row_perc" ] <- max(rowMeans(is.na(df_allPS)))
    df_summary[source_name,"PS_med_nan_row_perc" ] <- median(rowMeans(is.na(df_allPS)))

    # Plot hist of nan per junction in given data souce 
    ggplot(NULL, aes(x=rowMeans(is.na(df_allPS)))) + geom_histogram() + 
            scale_x_continuous(limits = c( 0,1 )) + xlab("Percent NAN in junction") +
            ggtitle(source_name) + geom_vline(xintercept = 0.5)
    ggsave(file.path(opt$out_dir,
            paste(out_prefix, source_name,"junction_nans", "png", sep = '.')),
            device = "png",
            width = 12,
            dpi = 300)
      }
  return(df_summary)
  }

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
      "source"=unique(df_metadata$data_source),
      "meta_file" = row[3]))
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
ls_source_names <- ls_metadata_files <- c()
for (item in ls_mesa_meta) {
     ls_mesa_inc_count_files <- append(ls_mesa_inc_count_files, item[1])
     ls_mesa_allPS_files <- append(ls_mesa_allPS_files, item[2])
     ls_mesa_cluster_files <- append(ls_mesa_cluster_files, item[3])
     ls_mesa_IR_table_files <- append(ls_mesa_IR_table_files, item[4])
     ls_mesa_IR_cov_dir <- append(ls_mesa_IR_cov_dir, item[5])
     ls_meta <- append(ls_meta, item[6])
     ls_sample_names <- append(ls_sample_names, item[7])
     ls_source_names <- append(ls_source_names, item[8]$source)
     ls_metadata_files <- append(ls_metadata_files, item[9])

   }

# Combine metadata from each data source by rows
df_merged_metadata <- do.call("rbind", ls_meta)
rownames(df_merged_metadata) <- c()

###############################
# Before nan adjust/filter
###############################

# # Run qc function on PS tables 
# df_summary_PS <- mesatable_qc(ls_mesa_allPS_files,  "PS")
# print(df_summary_PS)

# # Run qc function on IR tables
# df_summary_IR <- mesatable_qc(ls_mesa_IR_table_files,  "IR")
# print(df_summary_IR)

# For each df keep rows with less 20% nans then see how many overlapping nans remain 

#UMAP function labeling samples with high QC
make_umap_qc <- function( num_neighbor,meta_col, df, tag , df_meta, source ) {

  print(dim(df))



  set.seed(123)
  column.nan.counts <- colSums(is.na(df[,2:ncol(df)]))
  high_nan_samples = names(column.nan.counts[column.nan.counts > (dim(df)[1]/2)])

  #PCA
  prcomp.out = as.data.frame(prcomp(as.data.frame(t(df)), scale = F)$x)
  prcomp.out$Run = rownames(prcomp.out)
  prcomp.out.merge = merge(prcomp.out, y = df_meta)

  print(dim(df_meta))
  print(dim(prcomp.out))
  print(dim(prcomp.out.merge))


  #make palette
  n <- length(unique(df_meta[[meta_col]]))
  # print(n)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  #umap
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  print(dim(umap.out))

  #merge
  umap.out.merge = merge(umap.out, df_meta)

  print(dim(umap.out.merge))

  # print(head(umap.out.merge))
  print("gg:")


  #plot
  ggplot(umap.out.merge, aes(x, y, color = as.character(get(meta_col)))) +
          geom_point(size = 2) +
          theme_classic() +
          theme(legend.position="bottom",legend.title = element_blank()) +
          scale_color_manual(values=pal) +
          labs(title= paste("UMAP MESA: Cell types, n_neighbors =",num_neighbor,tag, sep = ' '))
          # +
          # geom_text(aes(x, y, label = "NAN"), data = umap.out.merge[umap.out.merge$Run %in% high_nan_samples,],vjust = 0, nudge_y = 0.25)

  #save
  ggsave(file.path(opt$out_dir,
          paste(source, "MESA_PCA_UMAP",meta_col,num_neighbor,tag,"nanqc.png", sep = '.')),
          device = "png",
          width = 12,
          dpi = 300)

  }

for (i in seq_along(ls_mesa_allPS_files)){
    source_name = ls_source_names[i]
    metadata = read.csv(file = as.character(ls_metadata_files[i]))

    print("------------------------------------")
    print(source_name)

    file_name = as.character(ls_mesa_allPS_files[i])
    if (!file.exists(file_name)){
      next
    }

    df_allPS <- read.table(as.character(file_name), row.names = 1, header=T)
    print("dim before nan drop:")
    print(dim(df_allPS))


    #keep junctions where less than 20% of samples are nans
    cutoff = ncol(df_allPS)*.20
    df_allPS_clean = df_allPS[rowSums(is.na(df_allPS))<cutoff,]
    print("dim after nan junction drop:")
    print(dim(df_allPS_clean))
    print("% na")
    print(sum(is.na(df_allPS_clean))/prod(dim(df_allPS)))

    #keep samples with less than 50% nan
    s_cutoff = nrow(df_allPS_clean)*.50
    print("dim after nan sample drop:")
    #keep col if number of nans is less than %nan
    df_allPS_clean = df_allPS_clean[,colSums(is.na(df_allPS_clean))<s_cutoff]
    print(dim(df_allPS_clean))
    print("% na after both drops")
    print(sum(is.na(df_allPS_clean))/prod(dim(df_allPS)))


    # replace remaining nan with row(junction) PS median
    indx <- which(is.na(df_allPS_clean), arr.ind = TRUE)
    df_allPS_clean[indx] <- apply(df_allPS_clean, 1, median, na.rm = TRUE)[indx[,"row"]]

    print("dim after nan replace with median:")
    print(dim(df_allPS_clean))

    print("nan % before ")
    print(sum(is.na(df_allPS))/prod(dim(df_allPS)))
    print("nan % after ")
    print(sum(is.na(df_allPS_clean))/prod(dim(df_allPS)))

    # print("make umap")

    if("group_label_qc" %in% colnames(metadata))
    {
      # lapply(c(5, 15, 25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata, source = source_name)
      # lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata)
      lapply(c(3, 5, 15, 25), make_umap_qc, meta_col="group_label_qc", df = df_allPS_clean, tag="nan_adjusted",  df_meta = metadata, source = source_name)
    }

     if("main_label_qc" %in% colnames(metadata))
    {
      # lapply(c(5, 15, 25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata, source = source_name)
      # lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata)
      lapply(c(3, 5, 15, 25), make_umap_qc, meta_col="main_label_qc", df = df_allPS_clean, tag="nan_adjusted",  df_meta = metadata, source = source_name)
    }
    

     if("donor_qc" %in% colnames(metadata))
    {
      # lapply(c(5, 15, 25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata, source = source_name)
      # lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata)
      lapply(c(3, 5, 15, 25), make_umap_qc, meta_col="donor_qc", df = df_allPS_clean, tag="nan_adjusted",  df_meta = metadata, source = source_name)
    }


     if("treatment_qc" %in% colnames(metadata))
    {
      # lapply(c(5, 15, 25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata, source = source_name)
      # lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_general", df = na.omit(df_allPS), tag="nan_dropped", df_meta = metadata)
      lapply(c(3, 5, 15, 25), make_umap_qc, meta_col="treatment_qc", df = df_allPS_clean, tag="nan_adjusted",  df_meta = metadata, source = source_name)
    }

}