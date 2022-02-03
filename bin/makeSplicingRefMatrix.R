#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)

# This script ...
# (1) Reads and plots the results from runMESAcompare.R , which does 1 vs all comparisons
#     using uses LM22 and LM6 cell sub types
# (2) Reads and plots the results from compareWithinType.R, which does 1 vs all comparisons
#     within T-cells, Monocytes/macrophages, etc


##########################
# Functions
##########################
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  #' Function to save pheatmaps to a pdf file
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

import_mesa_css <- function(filename, topN, plot_out_dir, css_dir, meta_col){
  #' Import results from MESA compare sample set script to get the top N
  #' significant events into lists
  #' @param filename -
  #' @param topN -
  #' @param plot_out_dir - path to output directory
  #' @param path tom esa compare sample set output directory
  #' @return  list of top N positive evnets, top negative events, top N
  #' negative and top N positive

  # Filename to string
  LM22_type <-  substr(filename, 1, nchar(filename)-4)
  print(LM22_type)

  # Make output directories
  if (!dir.exists(paste0(plot_out_dir,LM22_type))){
    dir.create(paste0(plot_out_dir,LM22_type),
     recursive = TRUE, showWarnings = TRUE)
  }

  # Open mesa css
  df <- read.table(
            file = paste0(css_dir, filename),
            sep="\t", header = TRUE)

  # Get top events by pvalue
  top_sig_by_pval <- df %>%
    dplyr::arrange(p.value)

  # Get top negative delta events
  top_sig_by_pval_negdelta <- top_sig_by_pval %>%
    dplyr::filter(delta < -.1 ) %>%
    dplyr::arrange(delta) %>%
    head(topN)%>%
    pull(event)

  print("top_sig_by_pval_negdelta:")
  print(top_sig_by_pval_negdelta)

  # Make plots for top negative events
  lapply(top_sig_by_pval_negdelta,  plot_event, cell_type = LM22_type,
        LM_type=meta_col, out_dir = plot_out_dir)

  # Get top positive delta events
  top_sig_by_pval_posdelta <- top_sig_by_pval %>%
    dplyr::filter(delta > .1 ) %>%
    dplyr::arrange(desc(delta)) %>%
    head(topN)%>%
    pull(event)

  # Make plots for top positive events
  lapply(top_sig_by_pval_posdelta,  plot_event, cell_type = LM22_type,
          LM_type=meta_col,
          out_dir = plot_out_dir )

  return(list(
              "top_pos" = top_sig_by_pval_posdelta,
              "top_neg" = top_sig_by_pval_negdelta,
              "top_neg_and_pos" = unlist(list(top_sig_by_pval_negdelta,top_sig_by_pval_posdelta)))
              )
}

plot_event <- function(sig_event, cell_type, out_dir, LM_type){
  #' Make jitter plot for a given event in all samples
  #'
  #' @param sig_event - string for the significant event to ploit
  #' @param cell_type - string of the cell type the event was significant in
  #' @param out_dir - path to output directory
  #' @return NA

  df <- all_PS_meta %>%
    tibble::rownames_to_column(var = "event")%>%
    dplyr::filter(event %in% list(paste0(sig_event), paste0(LM_type)))%>%
    t() %>%
    as.data.frame()

  df_ <- df # copy df
  colnames(df_) <- c( "PS", paste0(LM_type)) #add column names from first row

  df_ <- df_[-1,] %>% # drop first row
          dplyr::filter(PS != "NaN")  #drop samples with Nan

  df_$PS <- as.numeric(df_$PS)

  p <- ggplot( df_, aes(x = get(LM_type), y = PS, color=get(LM_type))) +
      # geom_violin() +
      geom_jitter(position=position_jitter(0.2), alpha = 0.5) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(legend.position = "None") +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "cell type")

  ggsave(plot = p, filename = paste0(out_dir,cell_type,"/",sig_event,".png"))

}

unpack_import_css_res <- function(ls_res){
  #' Unpack the results from the import_mesa_css() into lists
  #'
  #' @param ls_res output from
  #' @param list of top positive evnet, top negative events, top negative and positive

  ls_top_pos <- ls_top_neg <- ls_top_neg_and_pos  <- c()
  for (item in ls_res) {
       ls_top_pos <- append(ls_top_pos, item[1])
       ls_top_neg <- append(ls_top_neg, item[2])
       ls_top_neg_and_pos <- append(ls_top_neg_and_pos, item[3])
     }

  ls_top_pos<- unlist(ls_top_pos)
  ls_top_neg<- unlist(ls_top_neg)
  ls_top_neg_and_pos<- unlist(ls_top_neg_and_pos)

  return(list(ls_top_pos,ls_top_neg,ls_top_neg_and_pos))
  }


make_pheatmap <- function(ls_events, label, df_meta, df_PS){
  #' Make heatmap using the pheatmap package using the given events
  #' and data
  #'
  #' @param ls_events - list of events of interest to use in heatmap
  #' @param label - label to use in output file path
  #' @param df_meta - df of metadata with Run, val, and data_source, LM6, LM22 columns
  #' @param df_PS - df of MESA all PS file

  # Filter MESA all PS file to events of interest
  df_all_PS_sig_events <- df_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_events[[3]]) %>%
    dplyr::select(noquote(order(colnames(.)))) %>%
    tibble::column_to_rownames('event')

  for (val in list("LM22", "LM6")){
  # DF to label samples(columns) with labels
  df_sample_annotations <- df_meta %>%
    dplyr::select(Run, val, data_source) %>%
    dplyr::arrange(Run) %>%
    tibble::column_to_rownames("Run")

  stopifnot(rownames(df_sample_annotations) == colnames(df_all_PS_sig_events))

  heatmap_res <- pheatmap(
    main = paste0(" "),
    df_all_PS_sig_events,
    # scale = "row",
    show_rownames=F,
    show_colnames=F,
    na_col = "grey",
    annotation_col = df_sample_annotations)

  save_pheatmap_pdf(
    heatmap_res,
    paste0(opt$out_dir,"/ref_matrix/",label,"_",val,".pdf"))
}
}

import_mesa_to_heatmap<- function(ls_cell_types, top_n,  label, css_dir, meta_col){
  #' Import results from MESA compare_sample_sets runs within a broader cell type
  #' (e.g. within T-cells). Find the top signficant events, make event level
  #' plots, make heatmaps of the events. This function calls import_mesa_css(),
  #' unpack_import_css_res(), and make_pheatmap()
  #'
  #' @param ls_cell_types - list of LM22 cell-types that were compared in
  #' compareWithinType.R script
  #' @param top_n - integer; how many of the top splicing events to use
  #' @param label - string to use to represent cell type in output files
  #' @return ls_top_events - list containing 3 list - top positive events, top
  #' negative events, and top negative and positive


  # Get output files from compareWithinType script
  ls_css_file_names <- list.files(css_dir,pattern = ".tsv")

  print(ls_css_file_names)

  # For input cell type list , convert to filename
  ls_cell_types_file <- c()
  for (val in ls_cell_types){
    new_val <- paste0(gsub(" ","_", val), ".tsv")
    ls_cell_types_file <- append(ls_cell_types_file, new_val)
  }

  # Intersect with the files that exist (Not all will have a mesa css output )
  ls_css_file_names_cell_type  <- intersect(ls_css_file_names, ls_cell_types_file)

  print(ls_css_file_names_cell_type)

  #Import files, find top signficant events, plot each event
  ls_res <- lapply(
                          ls_css_file_names_cell_type,
                          topN=top_n,
                          import_mesa_css, meta_col =meta_col,
                          plot_out_dir =  paste0(opt$out_dir,"/ref_matrix/within_type/"),
                          css_dir =  css_dir)

  # Unpack top events into lists
  ls_top_events <- unpack_import_css_res(ls_res)

  # Filter metadata
  df_metadata_subset <- metadata %>%
    dplyr::filter(LM22 %in% ls_cell_types) %>%
    droplevels(.) %>%
    dplyr::arrange(Run)

  # Filter all PS
  df_all_PS <- all_PS %>%
    dplyr::select(as.vector(unlist(df_metadata_subset$Run)))

  # Make heatmap with this cell types events only within this cell types samples
  make_pheatmap(ls_top_events, paste0(label, "_diff_splicing_heatmap"),
          df_metadata_subset, df_all_PS )

  # Make heatmap with this cell types events and all samples
  make_pheatmap(ls_top_events, paste0(label,"_diff_splicing_heatmap_all_samples"),
          metadata, all_PS )

  return(ls_top_events)
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
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general) %>%
  tibble::column_to_rownames("Run") %>%
  t()

print(head(df_sample_annotations))

all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
print(head(all_PS))
print(dim(all_PS))

all_PS_meta <- rbind(all_PS, df_sample_annotations)
print(head(all_PS_meta))
print(tail(all_PS_meta))

print(dim(all_PS_meta))


# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/LM22"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/LM22"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/LM6"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/LM6"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/within_type"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/within_type"),
   recursive = TRUE, showWarnings = TRUE)
}

#########################################
# Import LM22 1 vs all comparisons
########################################
# Get all outputs from compare sample sets 1 vs all comparisons
ls_lm22_css_file_names <- list.files(
                        paste0(opt$out_dir,
                        "/compare_LM22/mesa_css_outputs/"),
                      pattern = ".tsv")
ls_lm22_css_file_names <- ls_lm22_css_file_names[!ls_lm22_css_file_names %in% c("heatmaps")]

# Import, find signficant events and plot each one
ls_lm22_res <- lapply(ls_lm22_css_file_names, topN=10,  import_mesa_css,
              meta_col  = "LM22",
                  plot_out_dir = paste0(opt$out_dir,"/ref_matrix/LM22/"),
                css_dir=paste0(
                  opt$out_dir,
                  "/compare_LM22/mesa_css_outputs/"))

print(ls_lm22_res)

# Unpack top events into lists
ls_lm22_top_events <- unpack_import_css_res(ls_lm22_res)
print(ls_lm22_top_events)

# Make heatmap using top events
make_pheatmap(ls_lm22_top_events, "LM22_diff_splicing_heatmap", metadata, all_PS )

#########################################
# Import LM6 1 vs all comparisons
########################################
# Get all outputs from compare sample sets 1 vs all comparisons
ls_lm6_css_file_names <- list.files(
                        paste0(opt$out_dir,
                        "/compare_LM6/mesa_css_outputs/"),
                      pattern = ".tsv")
ls_lm6_css_file_names <- ls_lm6_css_file_names[!ls_lm6_css_file_names %in% c("heatmaps")]

# Import, find signficant events and plot each one
ls_lm6_res <- lapply(ls_lm6_css_file_names, topN=10,  import_mesa_css,
              meta_col  = "LM6",
                  plot_out_dir = paste0(opt$out_dir,"/ref_matrix/LM6/"),
                css_dir=paste0(
                  opt$out_dir,
                  "/compare_LM6/mesa_css_outputs/"))

print(ls_lm6_res)

# Unpack top events into lists
ls_lm6_top_events <- unpack_import_css_res(ls_lm6_res)
print(ls_lm6_top_events)

# Make heatmap using top events
make_pheatmap(ls_lm6_top_events, "LM6_diff_splicing_heatmap", metadata, all_PS )

#########################################
# T-cell 1 vs all comparisons
########################################
# print("T-cells .........................")
# Get samples with this cell type
T_cell_types <- list(
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory  activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta")

# Import files, find top signficant events, make event level plots event, make heatmaps
ls_Tcell_top_events <- import_mesa_to_heatmap(
                T_cell_types, top_n=10, label= "Tcells",
                css_dir=paste0(
                  opt$out_dir,
                  "/compare_within_type/mesa_css_outputs/"),
                meta_col="LM22")
print(ls_Tcell_top_events)


######################################################
# Monocytes and macrophages 1 vs all comparisons
######################################################
print("Monocytes and macrophages .........................")

# Get samples with this cell type
mon_mac_cell_types <- list(
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2")

# Import files, find top signficant events, make event level plots event, make heatmaps
ls_mon_mac_top_events <- import_mesa_to_heatmap(
                          mon_mac_cell_types, top_n = 10,
                          label = "Monocytes_macrophages",
                          css_dir=paste0(opt$out_dir, "/compare_within_type/mesa_css_outputs/"),
                        meta_col="LM22")

print(ls_mon_mac_top_events)


##########################
# B-cells
##########################
# Get samples with this cell type
B_cell_types <- list(
  "B cells naive",
  "B cells memory")

# Import files, find top signficant events, make event level plots event, make heatmaps
ls_Bcell_top_events <- import_mesa_to_heatmap(
                          B_cell_types, top_n = 10,
                          label = "Bcells",
                          css_dir=paste0(opt$out_dir, "/compare_within_type/mesa_css_outputs/"),
                        meta_col="LM22")


##########################
# Dendritic cells
##########################
dendritic_cell_types <- list(
  "Dendritic cells resting",
  "Dendritic cells activated")

# Import files, find top signficant events, make event level plots event, make heatmaps
ls_dendritic_top_events <- import_mesa_to_heatmap(
                          dendritic_cell_types, top_n = 10,
                          label = "Dendritic",
                          css_dir=paste0(opt$out_dir, "/compare_within_type/mesa_css_outputs/"),
                        meta_col="LM22")


##########################
# Mast cells
##########################
# mast_cell_types <- list(
#   "Mast cells resting",
#   "Mast cells activated")
#
# # Import files, find top signficant events, make event level plots event, make heatmaps
# ls_mast_top_events <- import_mesa_to_heatmap(
#                           mast_cell_types, top_n = 10,
#                           label = "Mast",
#                           css_dir=paste0(opt$out_dir, "/compare_within_type/mesa_css_outputs/"))


# ##########################
# # NK cells
# ##########################
# NK_cell_types <- list(
#   "NK cells resting",
#   "NK cells activated")
#

# Import files, find top signficant events, make event level plots event, make heatmaps
# ls_NK_top_events <- import_mesa_to_heatmap(
#                           NK_cell_types, top_n = 10,
#                           label = "NK",
#                           css_dir=paste0(opt$out_dir, "/compare_within_type/mesa_css_outputs/"))
#
