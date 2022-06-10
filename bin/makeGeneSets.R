#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
# library(uwot)
# library(RColorBrewer)
library(foreach)
# library(doParallel)
# cl <- makeCluster(detectCores() - 1, outfile = "")
# registerDoParallel(cl)

import_deseq2 <- function(filename, topN, plot_out_dir, deseq2_dir, meta_col, comparison_label){
  #' Import results from MESA compare sample set script to get the top N
  #' significant events into lists
  #' @param filename -
  #' @param topN -
  #' @param plot_out_dir - path to output directory
  #' @param path tom esa compare sample set output directory
  #' @return  list of top N positive evnets, top negative events, top N
  #' negative and top N positive

  # Filename to string
  label_type <- stringr::str_trim(substr(filename, 1, nchar(filename)-4))

  # Open deseq2
  df_res <- read.csv(file = paste0(deseq2_dir,filename), header = TRUE)

  # Get top negative delta events
  df_topN_DEG_down_reg <- df_res %>%
      dplyr::filter(padj <= .05 ) %>%
      dplyr::filter(log2FoldChange < 0 ) %>%
      dplyr::arrange(padj) %>%
      head(topN) 
    #   %>%
    #   select(X)

  # Get top negative delta events
  df_topN_DEG_up_reg <- df_res %>%
      dplyr::filter(padj <= .05 ) %>%
      dplyr::filter(log2FoldChange > 0 ) %>%
      dplyr::arrange(padj) %>%
      head(topN) 
    #   %>%
    #   select(X)

return(list(label_type,
            droplevels(df_topN_DEG_up_reg$X),
            droplevels(df_topN_DEG_down_reg$X),
            df_topN_DEG_up_reg,
            df_topN_DEG_down_reg ))
}

import_deseq2_within<- function(ls_cell_types, topN,  label, deseq2_dir, meta_col){
  #' Import results from MESA compare_sample_sets runs within a broader cell type
  #' (e.g. within T-cells). Find the top signficant events, make event level
  #' plots, make heatmaps of the events. This function calls import_deseq2(),
  #' unpack_import_css_res(), and make_pheatmap()
  #'
  #' @param ls_cell_types - list of  cell-types that were compared in
  #' compareWithinType.R script
  #' @param topN - integer; how many of the top splicing events to use
  #' @param label - string to use to represent cell type in output files
  #' @return ls_top_events - list containing 3 list - top positive events, top
  #' negative events, and top negative and positive


  # Get output files from compareWithinType script
  ls_css_file_names <- list.files(deseq2_dir,pattern = ".csv")
  print(ls_css_file_names)
  print(label)


  # For input cell type list , convert to filename
  ls_cell_types_file <- c()
  for (val in ls_cell_types){
    new_val <- paste0(gsub(" ","_", val), ".csv")
    ls_cell_types_file <- append(ls_cell_types_file, new_val)
  }

  # Intersect with the files that exist (Not all will have a mesa css output )
  ls_css_file_names_cell_type  <- intersect(ls_css_file_names, ls_cell_types_file)

  #Import files, find top signficant events, plot each event
  ls_res <- lapply(
                    ls_css_file_names_cell_type,
                    topN=topN,
                    import_deseq2,
                    meta_col =meta_col,
                    plot_out_dir =  paste0(opt$out_dir,"/ref_matrix/within_group/"),
                    deseq2_dir =  deseq2_dir,
                    comparison_label = "within_group")

#   # If only 2 CSS files, they should have identical resuls
#   # (comparing A vs B then B) , so only return the results of one
#   if (length(ls_css_file_names_cell_type) <= 2) {
#       #only keep first result
#       df_res <- ls_res[[1]]
#   } else {
#       df_res <- dplyr::bind_rows(ls_res)
#   }

#   ls_top_events <- df_res$X

#   # Filter metadata
#   df_metadata_subset <- df_metadata %>%
#     dplyr::filter(main_label %in% ls_cell_types) %>%
#     droplevels(.) %>%
#     dplyr::arrange(Run)

#   # Filter all PS
#   df_log2tpm_batch_corrrected_subset <- df_log2tpm_batch_corrrected %>%
#     dplyr::select(as.vector(unlist(df_metadata_subset$Run)))

#   # Make heatmap with this cell types events only within this cell types samples
#   make_pheatmap(ls_top_events, paste0(label, "_gene_heatmap"),
#           df_metadata_subset, df_log2tpm_batch_corrrected_subset )

#   # Make heatmap with this cell types events and ALL samples
#   make_pheatmap(ls_top_events, paste0(label,"_gene_heatmap_all_samples"),
#           df_metadata, df_log2tpm_batch_corrrected )

  return(ls_res)
}


# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--input_path"),
    type = "character",
    default = NULL,
    help = " `"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to metadata")

    )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
ls_out_paths <- list("/gmt" )
for (path in ls_out_paths){

  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}
}

# Open files
df_metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- df_metadata %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

print(df_sample_annotations)

########################################
# Import main_label 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_main_label_de_file_names <- list.files(
                                  paste0(opt$input_path,"/compare_main_label/deseq2_outputs/"),
                                  pattern = ".csv")

ls_cell_types_file <- list()
for (val in unique(list(df_metadata$main_label))){
    new_val <- paste0(gsub(" ","_", val), ".csv")
    ls_cell_types_file <- append(ls_cell_types_file, new_val)
  }

ls_de_file_names_cell_type  <- unlist(intersect(ls_main_label_de_file_names, ls_cell_types_file))

# Import, find signficant events and plot each one
ls_main_label_res <- foreach(i=ls_de_file_names_cell_type,
                      .packages=c('magrittr','dplyr','ggplot2')) %do% {
    import_deseq2(
      filename = i,
      topN = 200,
      meta_col="main_label",
      plot_out_dir = paste0(opt$input_path,"/ref_matrix/main_label/"),
      deseq2_dir=paste0(opt$input_path,"/compare_main_label/deseq2_outputs/"),
      comparison_label = "main_label"
      )
  }

print(length(ls_main_label_res))

#Check its existence
if (file.exists(paste0(opt$out_dir,"/gene_set.gmt"))) {
  #Delete file if it exists
  file.remove(paste0(opt$out_dir,"/gene_set.gmt"))
}

outfile = file(paste0(opt$out_dir,"/gene_set.gmt"), open = 'a') # open in “a”ppend mode

for (i in ls_main_label_res){
  name = i[1][1]
  print(name)
  # print(name, "\t", "desc" )

  ls_UP <- i[[2]]
  cat(paste0(name,"_main_UP"),"\t", "desc" ,"\t", paste0(ls_UP, sep = "\t"), '\n', file = outfile, sep = '', append = TRUE)

  ls_DN <- i[[3]]
  cat(paste0(name,"_main_DN"),"\t","desc","\t" ,paste0(ls_DN, sep = "\t"), '\n', file = outfile, sep = '', append = TRUE)

}

ls_df_main_list = list()
for (i in ls_main_label_res){
  name = i[1][1]
  df_UP <- i[[4]] %>%
    as.data.frame() %>%
      mutate(set = paste0(name,"_main_UP"))
  df_DN <- i[[5]] %>%
      as.data.frame() %>%
      mutate(set = paste0(name,"_main_DN"))
  df_UP_DN <-  dplyr::bind_rows(df_UP, df_DN)

  ls_df_main_list[[as.character(name)]] <- df_UP_DN 
}

df_main_label <- do.call("rbind", ls_df_main_list)
print(dim(df_main_label))


########################################
# Import group_label 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_group_label_de_file_names <- list.files(
                                  paste0(opt$input_path,"/compare_group_label/deseq2_outputs/"),
                                  pattern = ".csv")

ls_group_label_cell_types_file <- list()
for (val in unique(list(df_metadata$group_label))){
    new_val <- paste0(gsub(" ","_", val), ".csv")
    ls_group_label_cell_types_file <- append(ls_group_label_cell_types_file, new_val)
  }

ls_de_file_names_cell_type_group_label  <- unlist(intersect(ls_group_label_de_file_names, ls_group_label_cell_types_file))

# Import, find signficant events and plot each one
ls_group_label_res <- foreach(i=ls_de_file_names_cell_type_group_label,
                      .packages=c('magrittr','dplyr','ggplot2')) %do% {
    import_deseq2(
      filename = i,
      topN = 200,
      meta_col="group_label",
      plot_out_dir = paste0(opt$out_dir,"/ref_matrix/group_label/"),
      deseq2_dir=paste0(opt$input_path,"/compare_group_label/deseq2_outputs/"),
      comparison_label = "group_label"
      )
  }

for (i in ls_group_label_res){
  name = i[1][1]
  print(name)
  # print(name, "\t", "desc" )

  ls_UP <- i[[2]]
  cat(paste0(name,"_group_UP"),"\t", "desc","\t" ,paste0(ls_UP, sep = "\t"), '\n', file = outfile, sep = '', append = TRUE)

  ls_DN <- i[[3]]
  cat(paste0(name,"_group_DN"),"\t","desc","\t" ,paste0(ls_DN, sep = "\t"), '\n', file = outfile, sep = '', append = TRUE)

}

ls_df_group_list = list()
for (i in ls_group_label_res){
  name = i[1][1]
  df_UP <- i[[4]] %>%
    as.data.frame() %>%
      mutate(set = paste0(name,"_group_UP"))
  df_DN <- i[[5]] %>%
      as.data.frame() %>%
      mutate(set = paste0(name,"_group_DN"))
  df_UP_DN <-  dplyr::bind_rows(df_UP, df_DN)

  ls_df_group_list[[as.character(name)]] <- df_UP_DN 
}

df_group_label <- do.call("rbind", ls_df_group_list)

########################################
# Import within type comparisons
# #######################################
# Match main cell type to group labels 
ls_group_cell_types <- unlist(unique(df_metadata[["group_label"]]))
ls_main_cell_types <- unlist(unique(df_metadata[["main_label"]]))

group2main <- list()
for (i in ls_group_cell_types){
  print(i)

  m <- df_metadata %>%
      filter(group_label == i) 
  if (length(unique(m$main_label)) > 1) {
      group2main<- append(group2main,list(list(i, unique(m$main_label)) ))
  }
}

# Run 
ls_within_res <- foreach(i=group2main,
                        .packages=c('magrittr','dplyr','ggplot2','pheatmap')) %do% {

      import_deseq2_within(
        ls_cell_types = i[[2]],
        topN = 200,
        label = i[[1]],
        deseq2_dir = paste0(opt$input_path, "/compare_within_group/deseq2_outputs/"),
        meta_col = "main_label")

      }

print(ls_within_res)
for (ls in ls_within_res){
  for (i in ls){
    name = i[[1]][1]
    print(name)

    ls_UP <- i[[2]]
    cat(paste0(name,"_within_UP"),"\t","desc","\t" ,paste0(ls_UP, sep = "\t"), '\n', file = outfile, sep = '', append = TRUE)

    ls_DN <- i[[3]]
    cat(paste0(name,"_within_DN"),"\t","desc","\t" ,paste0(ls_DN, sep = "\t"), '\n', file = outfile, sep = '', append = TRUE)

  }
}

ls_df_within_list = list()
for (ls in ls_within_res){
  for (i in ls){
    name = i[[1]][1]
    print(name)

    df_UP <- i[[4]] %>%
      as.data.frame() %>%
        mutate(set = paste0(name,"_within_UP"))
    df_DN <- i[[5]] %>%
        as.data.frame() %>%
        mutate(set = paste0(name,"_within_DN"))
    df_UP_DN <-  dplyr::bind_rows(df_UP, df_DN)

    ls_df_within_list[[as.character(name)]] <- df_UP_DN 

  }
}
df_within_label <- do.call("rbind", ls_df_within_list)
print(df_within_label)
print(dim(df_within_label))

#################################
# Count occurances across sets 
#####################################

df_all_sets  <- do.call("rbind", list(df_main_label, df_group_label, df_within_label ))
print(dim(df_all_sets))

write.csv(df_all_sets, paste0(opt$out_dir,"/df_gene_sets.csv"))


print(head(df_all_sets))
print(unique(df_all_sets$set))

gene_2_sets <- df_all_sets %>% 
  group_by(X) %>% 
  summarize(context = list(set)) %>%
  mutate(sets = map_chr(context, toString)) %>%
  select(X,sets) 

print(gene_2_sets)
# count genes in matrix
df_countbygene <- df_all_sets %>% 
  count(X, sort = TRUE) %>%
  inner_join(gene_2_sets, by="X")
write.csv(df_countbygene, paste0(opt$out_dir,"/countbygene.csv"))
print(df_countbygene)

