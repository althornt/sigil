#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)

list2heatmap <- function(ls_events,heatmap_title,out_file_desc,meta_col_to_use ){
  #' Make a heatmap using events from the given list and all PS data
  #' @param ls_events - list of events to use in heatmap
  #' @param heatmap_title - Title of heatmap "Splicing - X "
  #' @param out_file_desc - string to include in file name
  #' @param meta_col_to_use - column from the metadata file to use for cell type
  #' @return NA

  # Filter MESA all PS file to events of interest
  df_all_PS_sig_events <- all_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_events) %>%
    tibble::column_to_rownames('event')

  # DF to label samples(columns) with general and more specific labels if they exist
  if("sigil_cell_type_treatment" %in% colnames(metadata)){
    df_sample_annotations <- metadata %>%
      dplyr::select(Run,sigil_general, sigil_cell_type_treatment) %>%
      tibble::column_to_rownames("Run")
  }
  else {
    df_sample_annotations <- metadata %>%
      dplyr::select(Run,sigil_general) %>%
      tibble::column_to_rownames("Run")
  }

  # Create and save heatmap
  if (length(ls_events) >= 2){

      heatmap_res <- pheatmap(
        main = paste0(" Splicing - ", heatmap_title),
        df_all_PS_sig_events,
        scale = "row",
        show_rownames=F,
        show_colnames=F,
        na_col = "grey",
        annotation_col = df_sample_annotations)

      save_pheatmap_pdf(
        heatmap_res,
        paste0(opt$out_dir,
              "/mesa_compare_outputs/mesa_css_outputs/heatmaps/",
              meta_col_to_use,"_", out_file_desc,
              "_diff_splicing_heatmap.pdf"))
    }
}

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

  # Make manifest 1 - given cell type
  df_m1_main_cell_type <- metadata %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = df_m1_main_cell_type,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/",
            paste0(str_cell_type_val),".tsv"))

  # Make manifest 2 - all other cell types
  df_m2_others <- metadata %>%
    dplyr::filter(get(meta_col_to_use) != cell_type_val) %>%
    dplyr::select(Run)

  write.table(x = df_m2_others,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/not_",
            paste0(str_cell_type_val),".tsv"))

  file_css_out <- paste0(opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",
                               str_cell_type_val,".tsv")

  # If enough samples, compare groups
  if ((nrow(df_m1_main_cell_type)>2) & (nrow(df_m2_others)>2)){

    print("running MESA compare...")

    # Run MESA compare_sample_sets command ; 2>&1 sends standard error standard output
    cmd <- paste0(
      "mesa compare_sample_sets --psiMESA ",opt$mesa_PS," -m1 ",
      opt$out_dir,"/mesa_compare_outputs/manifests/",str_cell_type_val,".tsv",
      " -m2 ",opt$out_dir, "/mesa_compare_outputs/manifests/not_",str_cell_type_val,".tsv  -o",
      opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
      opt$gtf, " 2>&1")
    system(cmd)
  }

  # Read in MESA compare_sample_sets results if they were  made
  # (mesa doesnt run if not enough samples)


  if (file.exists(file_css_out)){
    df_css_out = read.csv(file = file_css_out, sep="\t")

    # Filter MESA all PS file to significant events
    ls_sig_css_out <- df_css_out %>%
      dplyr::filter((corrected < .01) & (abs(delta) > .2)) %>%
      dplyr::pull(event)

    list2heatmap(ls_sig_css_out,cell_type_val,str_cell_type_val,meta_col_to_use)
    return(ls_sig_css_out)

  } else{
    # return empty vector , since no diff splicing
    return( vector())
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

# Open files
metadata = read.csv(file = opt$metadata)
all_PS = read.table(file = opt$mesa_PS, sep="\t", row.names = 1, header = TRUE)
# mesa_ps = read.csv(file = opt$mesa_PS, sep="\t", row.names = "cluster")

# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/mesa_compare_outputs/manifests/"))){
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/manifests/"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/heatmaps/"))){
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/heatmaps/"),
   recursive = TRUE, showWarnings = TRUE)
}

###################
# sigil_general
###################
if("sigil_general" %in% colnames(metadata)){
  ls_general_cell_types <- unique(metadata[["sigil_general"]])

  print(ls_general_cell_types)

  # Run MESA compare_sample_sets for each general subtype and make heatmap
  sigil_general_mesa_comp_res <- sapply(
    ls_general_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="sigil_general",
    USE.NAMES = TRUE)

  ls_combined_diff_splice_events <- unlist(sigil_general_mesa_comp_res)
  print(length(ls_combined_diff_splice_events))

  list2heatmap(ls_combined_diff_splice_events,
              "All significant events",
              "all_diff_splicing_sigil_general",
              "sigil_general")
}

######################################
# sigil_cell_type_treatment
######################################

if("sigil_cell_type_treatment" %in% colnames(metadata)){
  ls_sigil_cell_type_treatment_cell_types <- unique(metadata[["sigil_cell_type_treatment"]])

  # Run MESA compare_sample_sets for each general subtype and make heatmap
  sigil_treatment_mesa_comp_res <- sapply(
    ls_sigil_cell_type_treatment_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="sigil_cell_type_treatment",
    USE.NAMES = TRUE)

  ls_combined_diff_splice_events_treatment <- unlist(sigil_treatment_mesa_comp_res)
  print(length(ls_combined_diff_splice_events_treatment))

  list2heatmap(ls_combined_diff_splice_events_treatment,
              "All significant events from comparisons using treatment",
              "all_diff_splicing_sigil_cell_type_treatment",
              "sigil_cell_type_treatment")
}
