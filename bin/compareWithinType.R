#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)


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

    print("running MESA compare_sample_sets...")

    # Run MESA compare_sample_sets command ; 2>&1 sends standard error standard output
    cmd <- paste0(
      "mesa compare_sample_sets --psiMESA ",PS_path,
      " -m1 ",opt$out_dir,"/manifests/",str_cell_type_val,".tsv",
      " -m2 ",opt$out_dir, "/manifests/not_",str_cell_type_val,".tsv  -o ",
      opt$out_dir, "/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
      opt$gtf, " 2>&1")

    system(cmd)

  } else {
    print("Not running MESA compare_sample_sets because not at least 3 samples in each group...")
  }
}

call_run_css_cell_type <- function(ls_cell_types, label){
  #' This function calls runCompareSampleSets_1_vs_all to compare within
  #' the given list of cell types
  #' @params ls_cell_types - list of cell type names within category
  #' @params label - string to use iin naming outputs

  # Get samples with this cell type
  ls_samples_types <- metadata %>%
    dplyr::filter(LM22 %in% ls_cell_types)

  # Keep events with over 75% samples with data in EACH subset type
  # ls_events_keep <- list()
  n <- 0
  for (val in ls_cell_types){
    print(val)
    print(n)

    df_samples <- metadata  %>%
      dplyr::filter(LM22 %in% val)

    # Fitler PS df to samples of this type
    all_PS_nan_filt_sub <- all_PS_nan_filt %>%
      dplyr::select(as.vector(unlist(df_samples$Run)))

    # Keep events with atleast 75% of data not missing in this cell type
    all_PS_nan_filt_sub_nans <- all_PS_nan_filt_sub[which(
                              rowMeans(!is.na(all_PS_nan_filt_sub)) >= 0.75), ]
    # Add list to list
    # print(length(ls_events_keep))

    # ls_events_keep <- append(
    #                   ls_events_keep,
    #                   rownames(all_PS_nan_filt_sub_nans)
    #                   )

    print("this type")
    print(length(rownames(all_PS_nan_filt_sub_nans)))

    if (n >0){
      print("intersect")
      print("before add")
      print(length(ls_events_keep))
      ls_events_keep <- Reduce(
                            intersect,
                            list(ls_events_keep, rownames(all_PS_nan_filt_sub_nans)))
    } else{
      print("else")
      ls_events_keep <-  rownames(all_PS_nan_filt_sub_nans)
    }


    print("after add")
    print(length(ls_events_keep))

    n <- n+ 1

    }

    print("final:")
    print(length(ls_events_keep))


    # print("intersection")
    # # Get intersection of events in each list
    # ls_events_keep_intersect <- Reduce(intersect,  ls_events_keep)
    #
    # print(length(ls_events_keep_intersect))
    # print(ls_events_keep_intersect)


  # Keep events that pass NAN filtering in ALL cell types within this category

  df_all_PS_nan_filt_subset <- all_PS_nan_filt %>%
    tibble::rownames_to_column(., "Run") %>%
    dplyr::filter(Run %in% ls_events_keep) %>%
    tibble::column_to_rownames(., "Run") %>%
    dplyr::select(as.vector(ls_samples_types$Run))

  # Write to file to be used by MESA compare
  path_all_PS_filt_out <- paste0(opt$out_dir, "/celltype_subset_dfs/",
                "mesa_allPS_nan_filt_",label,".tsv")
  write.table(x = df_all_PS_nan_filt_subset, na="nan", row.names = TRUE,
              quote=FALSE, col.names=NA, sep = "\t",
              file = path_all_PS_filt_out)

  # Run MESA compare_sample_sets using new reduced df (nan and sample filtered)
  # For each cell type within this general cell type (Tcells, mac and monoc, etc)
  sapply(
    unlist(ls_cell_types),
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="LM22",
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
if (!dir.exists(paste0(opt$out_dir,"/manifests/"))){
  dir.create(paste0(opt$out_dir,"/manifests/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/mesa_css_outputs/"))){
  dir.create(paste0(opt$out_dir,"/mesa_css_outputs"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/celltype_subset_dfs/"))){
  dir.create(paste0(opt$out_dir,"/celltype_subset_dfs/"),
   recursive = TRUE, showWarnings = TRUE)
}

# Open files
metadata = read.csv(file = opt$metadata)
all_PS = read.table(file = opt$mesa_PS, sep="\t", row.names = 1, header = TRUE)

# Remove rows with more than 50% NA
all_PS_nan_filt <- all_PS[which(rowMeans(!is.na(all_PS)) > 0.5), ]
write.table(x = all_PS_nan_filt,na="nan", row.names = TRUE, quote=FALSE,
          col.names=NA, sep = "\t",
          file = paste0(opt$out_dir, "/batch_corr_mesa_allPS_LM22_nan_filt.tsv"))
print("Number of junctions removed for having over 50% samples with Nans:")
print(nrow(all_PS)- nrow(all_PS_nan_filt))

# print(head(all_PS_nan_filt))


ls_lm22_cell_types <- unique(metadata[["LM22"]])
print(ls_lm22_cell_types)

##########################
# T-cells
##########################
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

# Do 1 vs all comparsino within T cell cell types
call_run_css_cell_type(T_cell_types, "Tcell" )

######################################################
# Monocytes and macrophages
######################################################

# Get samples with this cell type
# mon_mac_cell_types <- list(
#   "Monocytes",
#   "Macrophages M0",
#   "Macrophages M1",
#   "Macrophages M2")

# dont include ones without samples or else nan filtering will break

mon_mac_cell_types <- list(
    "Monocytes",
    "Macrophages M0",
    "Macrophages M1")

# Do 1 vs all comparsino within cell types
call_run_css_cell_type(mon_mac_cell_types, "Mon_Mac" )

##########################
# B-cells
##########################
# Get samples with this cell type
B_cell_types <- list(
  "B cells naive",
  "B cells memory")

# Do 1 vs all comparsino within cell types
call_run_css_cell_type(B_cell_types, "Bcell" )

##########################
# Dendritic cells
##########################
dendritic_cell_types <- list(
  "Dendritic cells resting",
  "Dendritic cells activated")

# Do 1 vs all comparsino within cell types
call_run_css_cell_type(dendritic_cell_types, "Dendritic" )

# ##########################
# # Mast cells
# ##########################
# mast_cell_types <- list(
#   "Mast cells resting",
#   "Mast cells activated")
#
# # Do 1 vs all comparsino within cell types
# call_run_css_cell_type(mast_cell_types, "Mast" )
#
# ##########################
# # NK cells
# ##########################
# NK_cell_types <- list(
#   "NK cells resting",
#   "NK cells activated")
#
# # Do 1 vs all comparsino within cell types
# call_run_css_cell_type(NK_cell_types, "NK" )
