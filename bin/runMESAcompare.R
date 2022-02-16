#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)


runCompareSampleSets_1_vs_all <- function(meta_col_to_use, cell_type_val){
  #' Run MESA compare_sample_sets comparing the given cell_type_val and
  #' create a heatmap of sifnificant events
  #'
  #' @param meta_col_to_use - column from the metadata file to use for cell type
  #' usually either "sigil_general" or "sigil_cell_type"
  #' @param cell_type_val - cell type name that is in the meta_col_to_use to
  #' compare to all other cell types
  #' @return List of differentially spliced events

  # Convert the given cell type to string with no spaces
  str_cell_type_val <- paste(unlist(strsplit(
                                    as.character(cell_type_val), split=" ")),
                                     collapse="_")

  print(str_cell_type_val)

  # Make manifest 1 - given cell type
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

  file_css_out <- paste0(opt$out_dir, "/compare_",meta_col_to_use,"/mesa_css_outputs/",
                               str_cell_type_val,".tsv")



  # If enough samples, compare groups
  if ((nrow(df_m1_main_cell_type)>2) & (nrow(df_m2_others)>2)){

    print("running MESA compare...")

    # Run MESA compare_sample_sets command ; 2>&1 sends standard error standard output
    # Can use batch_corr_mesa_allPS_LM22_nan_filt for both LM6 and LM22 comparisons
    cmd <- paste0(
      "mesa compare_sample_sets --psiMESA ",opt$out_dir, "/batch_corr_mesa_allPS_LM22_nan_filt.tsv",
      " -m1 ",opt$out_dir,"/compare_",meta_col_to_use,"/manifests/",str_cell_type_val,".tsv",
      " -m2 ",opt$out_dir, "/compare_",meta_col_to_use,"/manifests/not_",str_cell_type_val,".tsv  -o ",
      opt$out_dir, "/compare_",meta_col_to_use,"/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
      opt$gtf, " 2>&1")
    system(cmd)
    print(cmd)
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

if (!dir.exists(paste0(opt$out_dir,"/compare_LM6/manifests/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM6/manifests/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/compare_LM6/mesa_css_outputs/"))){
  dir.create(paste0(opt$out_dir,"/compare_LM6/mesa_css_outputs/"),
   recursive = TRUE, showWarnings = TRUE)
}

# Open files
metadata = read.csv(file = opt$metadata)
all_PS = read.table(file = opt$mesa_PS, sep="\t", row.names = 1, header = TRUE)

print("all_PS")
print(head(all_PS))
print(dim(all_PS))

# Remove rows with less than 25% NA
all_PS_nan_filt <- all_PS[which(rowMeans(!is.na(all_PS)) > 0.75), ]

print("all_PS_nan_filt")
print(dim(all_PS_nan_filt))

write.table(x = all_PS_nan_filt,na="nan", row.names = TRUE, quote=FALSE,
          col.names=NA, sep = "\t",
          file = paste0(opt$out_dir, "/batch_corr_mesa_allPS_LM22_nan_filt.tsv"))
print("Number of junctions removed for having over 75% samples with Nans:")
print(nrow(all_PS)- nrow(all_PS_nan_filt))

# print(sum(rowSums(is.na(all_PS))))

# print(head(metadata))
#
# bcells <- metadata %>%
#   dplyr::filter(LM6  == "B cells") %>%
#   droplevels()
#
# print(bcells$Run)
#
# print(all_PS %>%
# dplyr::select(as.vector(bcells$Run)) %>%
# dplyr:::filter(row.names(all_PS) %in% c("1:198696909-198699563:-"))
# )

###################
# LM22
###################
if("LM22" %in% colnames(metadata)){
  ls_lm22_cell_types <- unique(metadata[["LM22"]])

  print(ls_lm22_cell_types)

  # Run MESA compare_sample_sets for each general subtype and make heatmap
  sigil_lm22_mesa_comp_res <- sapply(
    ls_lm22_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="LM22",
    USE.NAMES = TRUE)

  ls_combined_diff_splice_events <- unlist(sigil_lm22_mesa_comp_res)
  print(length(ls_combined_diff_splice_events))

}

###################
# LM6
###################
if("LM6" %in% colnames(metadata)){
  ls_lm6_cell_types <- unique(metadata[["LM6"]])
  ls_lm6_cell_types <- ls_lm6_cell_types[ls_lm6_cell_types != ""]

  # Run MESA compare_sample_sets for each general subtype and make heatmap
  sigil_lm6_mesa_comp_res <- sapply(
    ls_lm6_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="LM6",
    USE.NAMES = TRUE)

  ls_lm6_combined_diff_splice_events <- unlist(sigil_lm6_mesa_comp_res)
  print(length(ls_lm6_combined_diff_splice_events))

}
