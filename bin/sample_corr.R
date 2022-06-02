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


# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "PS file"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to metadata"))

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Open files
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  arrange(Run) %>%
  dplyr::select(Run,main_label,group_label) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
# all_PS_meta <- rbind(all_PS, df_sample_annotations)

# Make output directories
ls_out_paths <- list("/outputs" )
for (path in ls_out_paths){

  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}
}


metadata <- metadata %>% 
    arrange(Run) %>%
    tibble::column_to_rownames("Run") %>%
    select(main_label, group_label, data_source)


print(metadata)


# Read in all MESA PS 
df_all_PS <- read.table(file = opt$i, sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric) 
df_all_PS <- df_all_PS[ ,order(names(df_all_PS)) ]

# Reduce

print(head(df_all_PS))


corr_heatmap <- function(df, tag){

  # Find sample correlation; ignoring nans
  df_corr <- cor(df, method = "spearman", use="complete.obs") %>% as.data.frame
  df_corr <- df_corr[ ,order(names(df_corr)) ]
  df_corr <- df_corr[order(row.names(df_corr)),]


  print(df_corr)

  # Make heatmap of sample correlations 

  ls_main_col = c(
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

  ls_group_col = c(

    "Monocytes macrophages" = "#808000",
    "Dendritic cells" = "#f58231",
    "T cells" = "#dcbeff",              
    "B cells" = "#469990",
    "NK cells" = "#ffe119",
     "Eosinophils" = "#42d4f4",
        "Neutrophils" = "#4363d8"

  )

  ha <- HeatmapAnnotation(
      df = metadata, 
      name = "main_label", 
      col = list("main_label"= ls_main_col, 
                  "group_label" = ls_group_col,
                  "data_source" = c("Choi" = "black",
                                  "Monaco" = "darkgrey" , 
                                  "Song" = "lightgrey")),
      show_legend = TRUE,
      which = "column",
      gp = gpar(col = NA)
      )

  ha_row <- rowAnnotation(
      df = metadata, 
      name = "main_label", 
      col = list("main_label"= ls_main_col, 
                  "group_label" = ls_group_col,
                "data_source" = c("Choi" = "black",
                                  "Monaco" = "darkgrey" , 
                                  "Song" = "lightgrey")),
      show_legend = TRUE,
      gp = gpar(col = NA)
      )

  ###############################################################################
  # clustered unscaled heatmap 
  png(file=paste0(opt$o,"/corr_heatmap_unscaled_clustered_", tag,".png"),
      width = 15,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht <- ComplexHeatmap::Heatmap(df_corr,
                                name="Spearman",
                                  show_row_names= TRUE,
                                  show_column_names = FALSE,
                                  row_names_gp = grid::gpar(fontsize =5),
                                  top_annotation=ha,
                                  right_annotation = ha_row,
                                  # heatmap_legend_param = list(legend_gp = gpar(fontsize = 3)),
                                  show_row_dend = FALSE,
                                  heatmap_legend_param = list(
                                  # legend_direction = "horizontal", 
                                  # legend_height = unit(1, "cm"),
                                  legend_gp = gpar(fontsize = 5))
                )

  draw(ht, merge_legend = TRUE)
  dev.off()
}

corr_heatmap(df_all_PS, "all_samples")

