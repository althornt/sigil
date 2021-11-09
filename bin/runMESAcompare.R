#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)


runCompareSampleSets_1_vs_all <- function(meta_col_to_use, cell_type_val){

  print(cell_type_val)

  # Convert cell type to string with no spaces
  str_cell_type_val <- paste(unlist(strsplit(as.character(cell_type_val), split=" ")), collapse="_")

  # Make manifest 1
  m1_main_cell_type <- metadata %>%
    dplyr::filter(get(meta_col_to_use) == cell_type_val) %>%
    dplyr::select(Run)

  print(dim(m1_main_cell_type))
  write.table(x = m1_main_cell_type,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/",
            paste0(str_cell_type_val),".tsv"))

  # Make manifest 2
  m2_others <- metadata %>%
    dplyr::filter(get(meta_col_to_use) != cell_type_val) %>%
    dplyr::select(Run)

  print(dim(m2_others))

  write.table(x = m2_others,row.names = FALSE, quote=FALSE,col.names=FALSE,
            file = paste0(opt$out_dir, "/mesa_compare_outputs/manifests/not_",
            paste0(str_cell_type_val),".tsv"))

  # Make MESA compare_sample_sets command ; 2>&1 sends standard error standard output
  cmd <- paste0(
    "mesa compare_sample_sets --psiMESA ",opt$mesa_PS," -m1 ",
    opt$out_dir,"/mesa_compare_outputs/manifests/",str_cell_type_val,".tsv",
    " -m2 ",opt$out_dir, "/mesa_compare_outputs/manifests/not_",str_cell_type_val,".tsv  -o",
    opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",str_cell_type_val,".tsv --annotation ",
    opt$gtf, " 2>&1")

  # Run MESA compare_sample_sets
  system(cmd)

  # Read in MESA compare_sample_sets results
  css_out_file <- paste0(opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/",
                     str_cell_type_val,".tsv")
  css_out = read.csv(file = css_out_file, sep="\t")

  # Filter to significant events
  sig_css_out <- css_out %>%
    dplyr::filter((corrected < .01) & (abs(delta) > .2))
  ls_sig_css_out <- sig_css_out %>%
    dplyr::pull(event)

  all_PS_sig_events <- all_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_sig_css_out) %>%
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

  if (length(ls_sig_css_out) >= 2){

      heatmap_res <- pheatmap(
        main = paste0(" Splicing ", cell_type_val),
        all_PS_sig_events,
        scale = "row",
        show_rownames=F,
        show_colnames=F,
        na_col = "grey",
        annotation_col = df_sample_annotations)

      save_pheatmap_pdf(
        heatmap_res,
        paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/heatmaps/",meta_col_to_use,"_", str_cell_type_val,"_diff_splicing_heatmap.pdf"))
    }

  return(ls_sig_css_out)

}


# Function to save pheatmaps to a file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
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
  print("make")
  dir.create(paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/heatmaps/"),
   recursive = TRUE, showWarnings = TRUE)
}

###################
# sigil_general
###################

ls_general_cell_types <- unique(metadata[["sigil_general"]])

if("sigil_general" %in% colnames(metadata)){

  # Run MESA compare_sample_sets for each general subtype and make heatmap
  sigil_general_mesa_comp_res <- sapply(
    ls_general_cell_types,
    runCompareSampleSets_1_vs_all,
    meta_col_to_use="sigil_general",
    USE.NAMES = TRUE)

  ls_combined_DEG <- unlist(sigil_general_mesa_comp_res)

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

  all_PS_sig_events <- all_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_combined_DEG) %>%
    tibble::column_to_rownames('event')

  if (length(all_PS_sig_events) >= 2){

      heatmap_res <- pheatmap(
        main = paste0(" All differentially spliced junctions "),
        all_PS_sig_events,
        scale = "row",
        show_rownames=F,
        show_colnames=F,
        na_col = "grey",
        annotation_col = df_sample_annotations)

      save_pheatmap_pdf(
        heatmap_res,
        paste0(opt$out_dir,"/mesa_compare_outputs/mesa_css_outputs/heatmaps/sigil_general_ALL_diff_splicing_heatmap.pdf"))
    }





}
