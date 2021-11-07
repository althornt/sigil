#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)


runCompareSampleSets_1_vs_all <- function(meta_col_to_use, cell_type_val){

  print(cell_type_val)

  str_cell_type_val <- paste(unlist(strsplit(as.character(cell_type_val), split=" ")), collapse="_")

  # Make manifest 1
  m1_main_cell_type <- metadata %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = m1_main_cell_type,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/",paste0(str_cell_type_val),".tsv"))

  # Make manifest 2
  m2_others <- metadata %>%
    dplyr::filter(get(meta_col_to_use) != cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = m2_others,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/not_",paste0(str_cell_type_val),".tsv"))

  # Make MESA compare_sample_sets command ; 2>&1 sends standard error standard output
  cmd <- paste0("mesa compare_sample_sets --psiMESA ",opt$mesa_PS," -m1 ",
                opt$out_dir,"/mesa_compare_outputs/manifests/",str_cell_type_val,".tsv",
                " -m2 ",opt$out_dir, "/mesa_compare_outputs/manifests/not_",str_cell_type_val,".tsv  -o",
                opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",str_cell_type_val,".tsv 2>&1")

  # Run MESA compare_sample_sets
  system(cmd)

  #system(noquote(cmd), wait=TRUE)

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
    help = "full path to put outputs"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Open files
metadata = read.csv(file = opt$metadata)
mesa_ps = read.csv(file = opt$mesa_PS, sep="\t", row.names = "cluster")

# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/mesa_compare_outputs/manifests/"))){
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/manifests/"), recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/"))){
  print("make")
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/"), recursive = TRUE, showWarnings = TRUE)
}

###################
# sigil_general
###################

ls_general_cell_types <- unique(metadata[["sigil_general"]])

if("sigil_general" %in% colnames(metadata)){
  sigil_general_mesa_comp_res <- sapply(
    ls_general_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="sigil_general")
}
