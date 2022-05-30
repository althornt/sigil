#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
# library(uwot)
library(RColorBrewer)
# library(foreach)
library(doParallel)
# cl <- makeCluster(detectCores() - 1, outfile = "")
# registerDoParallel(cl)
library(uwot)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(purrr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)


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
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
# all_PS_meta <- rbind(all_PS, df_sample_annotations)

# Make output directories
ls_out_paths <- list("/gmt" )
for (path in ls_out_paths){

  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}
}

gct_in <- "/mnt_/benchmark/all_data.gct"

dat.gct <- read.delim(file=gct_in, skip=2) %>%
  select(-Description) %>%
  column_to_rownames("Name") 

dat.gct <- t(scale(t(dat.gct))) # scale and center rows

metadata <- metadata %>% 
    tibble::column_to_rownames("Run") %>%
    select(main_label)

# https://sashamaps.net/docs/resources/20-colors/

ls_col = c(
    "B cells memory" = "#469990", #Teal                           
    "B cells naive"  = "#000075", #Navy

    "Eosinophils" = "#42d4f4",
    "Neutrophils" = "#4363d8",

    "Monocytes" = "#800000", #Black
    "Macrophages M0" = "#9A6324",           
    "Macrophages M1"  = "#808000", #Maroon

    "Dendritic cells resting" = "#f58231",
    "Dendritic cells activated" = "#e6194B",
                        
    "NK cells resting" = "#ffe119",

    "T cells gamma delta"   = "#fabed4", #pink
    "T cells follicular helper" = "#ffd8b1", #apricot
    "T cells regulatory"  = "#fffac8", #beige
    "T cells CD4 naive"   = "#aaffc3", #mint
    "T cells CD8"         = "#dcbeff" #lavender
)


png(file=paste0(opt$out_dir,"/heatmap.png"),
    width = 9,
    height    = 10,
    units     = "in",
    res       = 1200)

ha <- HeatmapAnnotation(
    df = metadata, 
    name = "main_label", 
    col = list("main_label"= ls_col),
    #  na_col = "grey",
    # annotation_legend_param = list(),
    show_legend = TRUE,
    which = "column",
    gp = gpar(col = NA)
    # border = FALSE,
    # gap = unit(1, "points"),

    # show_annotation_name = TRUE,
    # annotation_name_gp = gpar(),
    # annotation_name_offset = NULL,
    # annotation_name_side = ifelse(which == "column", "right", "bottom"),
    # annotation_name_rot = NULL,

    # annotation_height = NULL,
    # annotation_width = NULL,
    # height = NULL,
    # width = NULL,
    # simple_anno_size = ht_opt$simple_anno_size,
    # simple_anno_size_adjust = FALSE
    )

ht <- ComplexHeatmap::Heatmap(dat.gct,
                              name="Z-Score",
                                show_row_names= TRUE,
                                show_column_names = FALSE,
                                row_names_gp = grid::gpar(fontsize =5),
                                top_annotation=ha,
                                # heatmap_legend_param = list(legend_gp = gpar(fontsize = 3)),
                                show_row_dend = FALSE,
                                heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                # legend_height = unit(1, "cm"),
                                legend_gp = gpar(fontsize = 5))
               )

draw(ht)
dev.off()
