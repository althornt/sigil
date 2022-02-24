#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
# library(pheatmap)
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

runCompareSampleSets_1_vs_all <- function(meta_col_to_use, cell_type_val){
  #' Run MESA compare_sample_sets comparing the given cell_type_val and
  #' create a heatmap of sifnificant events
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

  #############################################################################
  # Filter out events where the main cell type has < 75% samples with data
  #############################################################################
  # Get samples with this given cell type
  df_samples <- metadata  %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val)

  all_PS_nan_filt_sub <- all_PS_nan_filt %>%
    dplyr::select(as.vector(unlist(df_samples$Run)))

  # Keep events with atleast 75% of data not missing in this cell type
  all_PS_nan_filt_sub_nans <- all_PS_nan_filt_sub[which(
                            rowMeans(!is.na(all_PS_nan_filt_sub)) >= 0.75), ]

  print("Number of events removed:")
  print(nrow(all_PS_nan_filt_sub)- nrow(all_PS_nan_filt_sub_nans))

  ls_events_keep <- rownames(all_PS_nan_filt_sub_nans)

  # Filter to PS file events that pass NAN filtering in main cell type
  df_all_PS_nan_filt_subset <- all_PS_nan_filt %>%
    tibble::rownames_to_column(., "Run") %>%
    dplyr::filter(Run %in% ls_events_keep) %>%
    tibble::column_to_rownames(., "Run") %>%
    tibble::rownames_to_column(., "cluster")

  # Write to file to be used as input by MESA compare
  path_all_PS_filt_out <- paste0(opt$out_dir,"/compare_",meta_col_to_use,
                                "/celltype_subset_dfs/mesa_allPS_nan_filt_",
                                str_cell_type_val,".tsv")

  write.table(x = df_all_PS_nan_filt_subset,
              na="nan", row.names = FALSE, quote=FALSE, sep = "\t",
              file = path_all_PS_filt_out)

  #############################################################################
  # Make sample manifests
  #############################################################################
  # Make manifest 1 - given main cell type
  df_m1_main_cell_type <- metadata %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = df_m1_main_cell_type,row.names = FALSE, quote=FALSE,
            file = paste0(opt$out_dir, "/compare_",meta_col_to_use,"/manifests/",
            paste0(str_cell_type_val),".tsv"))

  # Make manifest 2 - all other cell types
  df_m2_others <- metadata %>%
    dplyr::filter((get(meta_col_to_use) != cell_type_val) & (get(meta_col_to_use) != "") ) %>%
    dplyr::select(Run)

  write.table(x = df_m2_others,row.names = FALSE, quote=FALSE,
            file = paste0(opt$out_dir, "/compare_",meta_col_to_use,"/manifests/not_",
            paste0(str_cell_type_val),".tsv"))

  #############################################################################
  # Run compare sample sets
  #############################################################################
  # If enough samples, compare groups
  if ((nrow(df_m1_main_cell_type)>2) & (nrow(df_m2_others)>2)){

    print("running MESA compare...")

    # Run MESA compare_sample_sets command ; 2>&1 sends standard error standard output
    # Can use batch_corr_mesa_allPS_LM22_nan_filt for both LM6 and LM22 comparisons
    cmd <- paste0(
      "mesa compare_sample_sets --psiMESA ",path_all_PS_filt_out,
      " -m1 ",opt$out_dir,"/compare_",meta_col_to_use,"/manifests/",str_cell_type_val,".tsv",
      " -m2 ",opt$out_dir, "/compare_",meta_col_to_use,"/manifests/not_",str_cell_type_val,".tsv  -o ",
      opt$out_dir, "/compare_",meta_col_to_use,"/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
      opt$gtf, " 2>&1")

    system(cmd)
    # print(cmd)
  }

}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  #' Function to save pheatmaps to a pdf file
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
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
if (!dir.exists(paste0(opt$out_dir,"/compare_LM22/manifests/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM22/manifests/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/compare_LM22/mesa_css_outputs/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM22/mesa_css_outputs/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/compare_LM22/celltype_subset_dfs/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM22/celltype_subset_dfs/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/compare_LM6/manifests/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM6/manifests/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/compare_LM6/mesa_css_outputs/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM6/mesa_css_outputs/"),
   recursive = TRUE, showWarnings = TRUE)
}
if (!dir.exists(paste0(opt$out_dir,"/compare_LM6/celltype_subset_dfs/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM6/celltype_subset_dfs/"),
   recursive = TRUE, showWarnings = TRUE)
}
# Open files
metadata = read.csv(file = opt$metadata)
all_PS = read.table(file = opt$mesa_PS, sep="\t", row.names = 1, header = TRUE)

print("all_PS")
# print(head(all_PS))
print(dim(all_PS))

# Remove rows with more than 25% NA
all_PS_nan_filt <- all_PS[which(rowMeans(!is.na(all_PS)) > 0.75), ]

print("all_PS_nan_filt")
print(dim(all_PS_nan_filt))

# Write table and add name rownames "cluster"
write.table(x = data.frame("cluster"=rownames(all_PS_nan_filt),all_PS_nan_filt),
          na="nan", row.names = FALSE, quote=FALSE,
           sep = "\t",
          file = paste0(opt$out_dir, "/batch_corr_mesa_allPS_LM22_nan_filt.tsv"))
print("Number of junctions removed for having over 75% samples with Nans:")
print(nrow(all_PS)- nrow(all_PS_nan_filt))

# read_in_all_PS_nan_filt = read.table(file = paste0(
#       opt$out_dir, "/batch_corr_mesa_allPS_LM22_nan_filt.tsv"),
#        sep="\t", row.names = 1, header = TRUE)
#
#
# print("read_in_all_PS_nan_filt after read in ")
# print(head(read_in_all_PS_nan_filt))
# print(dim(read_in_all_PS_nan_filt))

print(cl)

###################
# LM22
###################
if("LM22" %in% colnames(metadata)){
  ls_lm22_cell_types <- unique(metadata[["LM22"]])

  foreach(i=ls_lm22_cell_types, .packages='magrittr') %dopar% {
    runCompareSampleSets_1_vs_all(
        cell_type_val = i,
        meta_col_to_use="LM22")}

}


##################
# LM6
##################
if("LM6" %in% colnames(metadata)){
  ls_lm6_cell_types <- unique(metadata[["LM6"]])
  ls_lm6_cell_types <- ls_lm6_cell_types[ls_lm6_cell_types != ""]

  foreach(i=ls_lm6_cell_types, .packages='magrittr') %dopar% {
    runCompareSampleSets_1_vs_all(
        cell_type_val = i,
        meta_col_to_use="LM6")}
}


###################
# LM22
###################
# if("LM22" %in% colnames(metadata)){
#   ls_lm22_cell_types <- unique(metadata[["LM22"]])
#
#   # Run MESA compare_sample_sets for each general subtype and make heatmap
#   sigil_lm22_mesa_comp_res <- sapply(
#     ls_lm22_cell_types,
#     runCompareSampleSets_1_vs_all,
#     meta_col_to_use="LM22",
#     USE.NAMES = TRUE)
#
#   ls_combined_diff_splice_events <- unlist(sigil_lm22_mesa_comp_res)
#   print(length(ls_combined_diff_splice_events))
#
# }

##################
# LM6
##################
# if("LM6" %in% colnames(metadata)){
#   ls_lm6_cell_types <- unique(metadata[["LM6"]])
#   ls_lm6_cell_types <- ls_lm6_cell_types[ls_lm6_cell_types != ""]
#
#   # Run MESA compare_sample_sets for each general subtype and make heatmap
#   sigil_lm6_mesa_comp_res <- sapply(
#     ls_lm6_cell_types,
#     runCompareSampleSets_1_vs_all,
#     meta_col_to_use="LM6",
#     USE.NAMES = TRUE)
#
#
# }
