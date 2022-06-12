#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(uwot)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(purrr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "directory to "),

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

# Make output directories
if (!dir.exists(paste0(opt$out_dir))){
  dir.create(paste0(opt$out_dir),
  recursive = TRUE, showWarnings = TRUE)
}


# Open files
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  arrange(Run) %>%
  dplyr::select(Run,main_label,group_label, sigil_general) %>%
  # dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Read in and sort metadata 
metadata <- metadata %>% 
    arrange(Run) %>%
    tibble::column_to_rownames("Run") %>%
    select(main_label, group_label)
    # select(main_label, group_label, data_source)

#############
# main label
############
ls_dir<- list.dirs(path = opt$i, full.names = FALSE, recursive = TRUE)

ls_main_dfs <- list()
for (d in ls_dir){

  filename <- paste0(opt$i,d,"/eval_by_main_label.csv")

  if (file.exists(filename)){
    df_enr <- read.csv(file = filename, stringsAsFactors = FALSE) %>%
      select(main_label, Mean) %>%
      rename(!!d := Mean) 
  
    ls_main_dfs <- append(df_enr, ls_main_dfs  )
  }
}

df_main_eval <- dplyr::bind_rows(ls_main_dfs, .id= "m") %>%
  as.data.frame() %>%
  select(-m) %>%
  tibble::column_to_rownames("main_label") 

df_main_eval[] <- sapply(df_main_eval, as.numeric)
  
png(file=paste0(opt$out_dir,"/heatmap_main.png"),
    width = 9,
    height    = 10,
    units     = "in",
    res       = 1200)

ht_order <- ComplexHeatmap::Heatmap(df_main_eval,
                              # name="Z-Score",
                                # column_order =order(colnames(df_main_eval)),
                                # row_order = order(rownames(df_main_eval)), 
                                show_row_names= TRUE,
                                show_column_names = TRUE,
                                row_names_gp = grid::gpar(fontsize =5),
                                column_names_gp = grid::gpar(fontsize =5),
                                # top_annotation=ha,
                                # heatmap_legend_param = list(legend_gp = gpar(fontsize = 3)),
                                show_row_dend = TRUE,
                                heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                # legend_height = unit(1, "cm"),
                                legend_gp = gpar(fontsize = 5))
              )

draw(ht_order, merge_legend = TRUE)
dev.off()

#############
# group label
############

ls_dir<- list.dirs(path = opt$i, full.names = FALSE, recursive = TRUE)

ls_group_dfs <- list()
df = data.frame(columns = unique(metadata$group_label))

for (d in ls_dir){
  filename <- paste0(opt$i,d,"/eval_by_group_label.csv")

  if (file.exists(filename)){
    df_enr <- read.csv(file = filename, stringsAsFactors = FALSE) %>%
      select(group_label, Mean) %>%
      rename(!!d := Mean) 
  
    ls_group_dfs <- append(df_enr, ls_group_dfs  )
  }
}

df_group_eval <- dplyr::bind_rows(ls_group_dfs, .id= "m") %>%
  as.data.frame() %>%
  select(-m) %>%
  tibble::column_to_rownames("group_label") 
df_group_eval[] <- sapply(df_group_eval, as.numeric)
  
png(file=paste0(opt$out_dir,"/heatmap_group.png"),
    width = 9,
    height    = 10,
    units     = "in",
    res       = 1200)

ht_order <- ComplexHeatmap::Heatmap(df_group_eval,
                                # col = colorRamp2(c(0, 1), c("white", "red")), 
                              # name="Z-Score",
                                # column_order =order(colnames(df_group_eval)),
                                # row_order = order(rownames(df_group_eval)), 
                                show_row_names= TRUE,
                                show_column_names = TRUE,
                                row_names_gp = grid::gpar(fontsize =5),
                                column_names_gp = grid::gpar(fontsize =5),
                                # top_annotation=ha,
                                # heatmap_legend_param = list(legend_gp = gpar(fontsize = 3)),
                                show_row_dend = TRUE,
                                heatmap_legend_param = list(
                                # legend_direction = "horizontal", 
                                # legend_height = unit(1, "cm"),
                                legend_gp = gpar(fontsize = 5))
              )

draw(ht_order, merge_legend = TRUE)
dev.off()