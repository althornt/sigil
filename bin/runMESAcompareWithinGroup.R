#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

#############################
# Functions
##########################
runCompareSampleSets_1_vs_all <- function(meta_col_to_use, ls_gen_cell_types,
                                                  cell_type_val, PS_path){
  #' Run MESA compare_sample_sets comparing the given cell_type_val
  #'
  #' @param meta_col_to_use - column from the metadata file to use for cell type
  #' usually either "sigil_general" or "sigil_cell_type"
  #' @param ls_gen_cell_types - list of cell type names within category
  #' @param cell_type_val - cell type name that is in the meta_col_to_use to
  #' compare to all other cell types
  #' @param PS_path - path to filtered all PS file to use

  # Convert the given cell type to string with no spaces
  str_cell_type_val <- paste(unlist(strsplit(
                                    as.character(cell_type_val), split=" ")),
                                    collapse="_")

  print(str_cell_type_val)

  # Make manifest 1 - given cell type
  df_m1_main_cell_type <- metadata %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = df_m1_main_cell_type,
            row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/manifests/",
            paste0(str_cell_type_val),".tsv"))

  # Make manifest 2 - all other cell types
  df_m2_others <- metadata %>%
    dplyr::filter(get(meta_col_to_use) %in% ls_gen_cell_types) %>%
    dplyr::filter(get(meta_col_to_use) != cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = df_m2_others,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/manifests/not_",
            paste0(str_cell_type_val),".tsv"))

  # If enough samples, compare groups
  if ((nrow(df_m1_main_cell_type)>2) & (nrow(df_m2_others)>2)){

    # Run MESA compare_sample_sets command ; 2>&1 sends standard error standard output
    cmd <- paste0(
      "mesa compare_sample_sets --psiMESA ",PS_path,
      " -m1 ",opt$out_dir,"/manifests/",str_cell_type_val,".tsv",
      " -m2 ",opt$out_dir, "/manifests/not_",str_cell_type_val,".tsv  -o ",
      opt$out_dir, "/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
      opt$gtf, " 2>&1")

    system(cmd)
    print(cmd)

  } else {
    print("Not running MESA compare_sample_sets because not at least 3 samples in each group...")
  }
}

call_run_css_cell_type <- function(ls_cell_types, label){
  #' This function calls runCompareSampleSets_1_vs_all to compare within
  #' the given list of cell types
  #' @param ls_cell_types - list of cell type names within category
  #' @param label - string to use iin naming outputs

  print(label)

   # Convert the given cell type to string with no spaces
  label <- paste(unlist(strsplit(
                                    as.character(label), split=" ")),
                                    collapse="_")

  ls_cell_types <- unlist(ls_cell_types)
  print(ls_cell_types )

  # Get samples with this given cell type
  ls_samples_types <- metadata %>%
    dplyr::filter(main_label %in% ls_cell_types)

  # Keep events with over 75% samples with data in EACH subset type
  # Loop through subsets to get the passing events and find intersect
  n <- 0
  for (val in ls_cell_types){
    # Filter metadata and PS df to samples of this type
    df_samples <- metadata  %>%
      dplyr::filter(main_label %in% val)
    all_PS_nan_filt_sub <- all_PS_nan_filt %>%
      dplyr::select(as.vector(unlist(df_samples$Run)))

    # Keep events with atleast 75% of data not missing in this cell type
    all_PS_nan_filt_sub_nans <- all_PS_nan_filt_sub[which(
                              rowMeans(!is.na(all_PS_nan_filt_sub)) >= 0.75), ]

    if (n >0){
      ls_events_keep <- Reduce(
                          intersect,
                          list(ls_events_keep, rownames(all_PS_nan_filt_sub_nans)))
    } else{
      ls_events_keep <-  rownames(all_PS_nan_filt_sub_nans)
    }
      n <- n+ 1
    }

  # Filter to events that pass NAN filtering in ALL cell types within this category
  df_all_PS_nan_filt_subset <- all_PS_nan_filt %>%
    tibble::rownames_to_column(., "Run") %>%
    dplyr::filter(Run %in% ls_events_keep) %>%
    tibble::column_to_rownames(., "Run") %>%
    dplyr::select(as.vector(ls_samples_types$Run)) %>%
    tibble::rownames_to_column(., "cluster")

  # Write to file to be used by MESA compare
  path_all_PS_filt_out <- paste0(opt$out_dir, "/celltype_subset_dfs/",
                "mesa_allPS_nan_filt_",label,".tsv")

  write.table(x = df_all_PS_nan_filt_subset,
              na="nan", row.names = FALSE, quote=FALSE, sep = "\t",
              file = path_all_PS_filt_out)

  # Run MESA compare_sample_sets using new reduced df (nan and sample filtered)
  # For each cell type within this general cell type (Tcells, mac and monoc, etc)
  sapply(
    unlist(ls_cell_types),
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="main_label",
    ls_gen_cell_types=unlist(ls_cell_types),
    PS_path = path_all_PS_filt_out,
    USE.NAMES = TRUE)

}

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--mesa_PS"),
    type = "character",
    default = NULL,
    help = " `mesa_allPS.tsv` input file "),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-g", "--gtf"),
    type = "character",
    default = NULL,
    help = "full path to gtf file used in MESA"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
ls_out_paths <- list(
          "/manifests/", 
          "/mesa_css_outputs/",
          "/celltype_subset_dfs/"
          )
for (path in ls_out_paths){
  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
    }
}

# Open files
metadata = read.csv(file = opt$metadata, stringsAsFactors=FALSE)
df_all_PS = read.table(file = opt$mesa_PS, sep="\t", row.names = 1, header = TRUE)
print(dim(df_all_PS))

# Write filtered PS df and add name rownames "cluster" to keep expected mesa format
str_PS_basename <- basename(opt$mesa_PS)
str_PS_basename<- substr(str_PS_basename,1,nchar(str_PS_basename)-4)

# Remove rows with more than 50% NA
all_PS_nan_filt <- df_all_PS[which(rowMeans(!is.na(df_all_PS)) > 0.5), ]

# Check PS and filtered PS dfs have expected number of samples
if(all.equal(length(metadata$Run),
            ncol(df_all_PS),
            ncol(all_PS_nan_filt)) != TRUE)
            stop("Error: Number of columns(samples) is not consistent")

write.table(x = all_PS_nan_filt, na="nan", row.names = TRUE, quote=FALSE,
          col.names=NA, sep = "\t",
          file = paste0(opt$out_dir,"/", str_PS_basename,"_nan_filt.tsv"))
print("Number of junctions removed for having over 50% samples with Nans:")
print(nrow(df_all_PS)- nrow(all_PS_nan_filt))


# Match main cell type to group labels 
ls_group_cell_types <- unlist(unique(metadata[["group_label"]]))
ls_main_cell_types <- unlist(unique(metadata[["main_label"]]))

group2main <- list()
for (i in ls_group_cell_types){
  print(i)
  m <- metadata %>%
      filter(group_label == i) 
  if (length(unique(m$main_label)) > 1) {
      group2main<- append(group2main,list(list(i, unique(m$main_label)) ))
  }
}

print(group2main)

# ####################################################
# # Run 1 vs all within each group
# ####################################################

# Run in parallel
foreach(i=group2main, .packages=c('magrittr','dplyr')
  )%dopar%{
          call_run_css_cell_type(
            ls_cell_types = i[2] , 
            label = i[1]
            )
      }
