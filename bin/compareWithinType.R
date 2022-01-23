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
runCompareSampleSets_1_vs_all <- function(meta_col_to_use, cell_type_val, PS_path){
  #' Run MESA compare_sample_sets comparing the given cell_type_val
  #'
  #' @param meta_col_to_use - column from the metadata file to use for cell type
  #' usually either "sigil_general" or "sigil_cell_type"
  #' @param cell_type_val - cell type name that is in the meta_col_to_use to
  #' compare to all other cell types

  # Convert the given cell type to string with no spaces
  str_cell_type_val <- paste(unlist(strsplit(
                                    as.character(cell_type_val), split=" ")),
                                     collapse="_")

  print(str_cell_type_val)

  # Make manifest 1 - given cell type
  df_m1_main_cell_type <- metadata %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val) %>%
    dplyr::select(Run)


  print(df_m1_main_cell_type)

  write.table(x = df_m1_main_cell_type,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/",
            paste0(str_cell_type_val),".tsv"))

  # Make manifest 2 - all other cell types
  df_m2_others <- metadata %>%
    dplyr::filter(get(meta_col_to_use) != cell_type_val) %>%
    dplyr::select(Run)

  print(df_m2_others)
  write.table(x = df_m2_others,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/not_",
            paste0(str_cell_type_val),".tsv"))

  # file_css_out <- paste0(opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",
  #                              str_cell_type_val,".tsv")

  # If enough samples, compare groups
  if ((nrow(df_m1_main_cell_type)>2) & (nrow(df_m2_others)>2)){

    print("running MESA compare...")

    # Run MESA compare_sample_sets command ; 2>&1 sends standard error standard output
    cmd <- paste0(
      "mesa compare_sample_sets --psiMESA ",PS_path,
      " -m1 ",opt$out_dir,"/mesa_compare_outputs/manifests/",str_cell_type_val,".tsv",
      " -m2 ",opt$out_dir, "/mesa_compare_outputs/manifests/not_",str_cell_type_val,".tsv  -o ",
      opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
      opt$gtf, " 2>&1")
    print(cmd)
    system(cmd)

  }
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
if (!dir.exists(paste0(opt$out_dir,"/mesa_compare_outputs/manifests/"))){
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/manifests/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/"))){
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs"),
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
          col.names=TRUE, sep = "\t",
          file = paste0(opt$out_dir, "/LM22_mesa_allPS_nan_filt.tsv"))
print("Number of junctions removed for having over 50% samples with Nans:")
print(nrow(all_PS)- nrow(all_PS_nan_filt))

ls_lm22_cell_types <- unique(metadata[["LM22"]])
print(ls_lm22_cell_types)

##########################
# T-cell
##########################
# T_cell_types <- list(
#   "T cells CD8",
#   "T cells CD4 naive",
#   "T cells CD4 memory resting",
#   "T cells CD4 memory  activated",
#   "T cells follicular helper",
#   "T cells regulatory (Tregs)",
#   "T cells gamma delta")
#
# print(metadata)
# # Get samples with this cell type
# ls_samples_T_cell_types <- metadata %>%
#   dplyr::filter(LM22 %in% T_cell_types) %>%
#   dplyr::pull(Run)
#
# # Reduce df to T-cell types ; write df to file
# Tcell_all_PS_nan_filt <- all_PS_nan_filt[ls_samples_T_cell_types]
#
# # Remove rows with more than 50% NA within Tcells
# Tcell_all_PS_nan_filt <- Tcell_all_PS_nan_filt[which(rowMeans(!is.na(all_PS)) > 0.5), ]
#
# #!!!! writing wrong file need to fix !!!
# write.table(x = all_PS_nan_filt, na="nan", row.names = TRUE, quote=FALSE,
#           col.names=TRUE, sep = "\t",
#           file = paste0(opt$out_dir, "/celltype_subset_dfs/","LM22_mesa_allPS_nan_filt_Tcells.tsv"))
#
#
# # Run MESA compare_sample_sets within each subtype using new reduced df
# sigil_lm22_mesa_comp_res <- sapply(
#   T_cell_types,
#   runCompareSampleSets_1_vs_all,
#   meta_col_to_use="LM22",
#   PS_path = paste0(opt$out_dir, "/celltype_subset_dfs/","LM22_mesa_allPS_nan_filt_Tcells.tsv"),
#   USE.NAMES = TRUE)

# ls_combined_diff_splice_events <- unlist(sigil_lm22_mesa_comp_res)
# print(length(ls_combined_diff_splice_events))

#
# ######################################################
# # Monocytes and macrophages
# ######################################################

# Get samples with this cell type
mon_mac_cell_types <- list(
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2")


call_run_css_cell_type <- function(ls_cell_types, label){
  #' This function calls runCompareSampleSets_1_vs_all to compare within
  #' the given list of cell types
  #' @params

  # Get samples with this cell type
  ls_samples_types <- metadata %>%
    dplyr::filter(LM22 %in% ls_cell_types) %>%
    dplyr::pull(Run)

  # Reduce df to T-cell types ; write df to file
  df_all_PS_nan_filt_subset <- all_PS_nan_filt[ls_samples_types]

  # Remove rows with more than 50% NA within Monocytes and Macrophages samples
  df_all_PS_nan_filt_subset <- df_all_PS_nan_filt_subset[which(rowMeans(!is.na(all_PS)) > 0.5), ]

  #!!!! writing wrong file need to fix !!!

  write.table(x = all_PS_nan_filt, na="nan", row.names = TRUE, quote=FALSE,
            col.names=TRUE, sep = "\t",
            file = paste0(opt$out_dir, "/celltype_subset_dfs/","mesa_allPS_nan_filt_",label,".tsv"))


  # Run MESA compare_sample_sets within each subtype using new reduced df
  sapply(
    ls_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="LM22",
    PS_path = paste0(opt$out_dir, "/celltype_subset_dfs/","mesa_allPS_nan_filt_",label,".tsv"),
    USE.NAMES = TRUE)


}

# Do 1 vs all comparsino within Monocytes and Macrophages cell types
call_run_css_cell_type(mon_mac_cell_types, "Mon_Mac" )
