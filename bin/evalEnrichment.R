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

ssgsea_acc <- function(dat.gct, out_dir, name){

  dat.gct_scale <- t(scale(t(dat.gct))) # scale and center rows

  # Drop gene sets that arent UP or DN
  dat.gct_UPDN <- dat.gct_scale %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Run") %>%
      dplyr::filter( grepl('UP|DN', Run)) %>%
      tibble::column_to_rownames("Run") 

  # Add row for which set has max and min score in each sample 
  max <- rownames(dat.gct_UPDN)[apply(dat.gct_UPDN,2,which.max)]
  dat.gct_UPDN["max",] <- max
  min <- rownames(dat.gct_UPDN)[apply(dat.gct_UPDN,2,which.min)]
  dat.gct_UPDN["min",] <- min

  dat.gct_UPDN["main_label", ] <- df_sample_annotations["main_label",]
  dat.gct_UPDN["data_source", ] <- df_sample_annotations["data_source",]

  # Transpose to compare label to max and min 
  # check if max value is label + "UP"
  # check if min value is label + "DN"
  df_ <- dat.gct_UPDN %>%
      t() %>%
      as.data.frame() %>%
      select(main_label,data_source, min, max) %>%
      mutate(main_label_str = gsub(" ", "_", main_label)) %>%
      mutate(max_match= ifelse((startsWith(as.character(max),
                                          as.character(main_label_str)) & 
                                  endsWith(as.character(max),
                                          as.character("UP"))) , 1, 0)) %>%
      mutate( min_match= ifelse((startsWith(as.character(min),
                                          as.character(main_label_str)) & 
                                  endsWith(as.character(min),
                                          as.character("DN"))) , 1, 0)) 
      
  # print(df_)

  print("min match")
  print(sum(df_$min_match)/nrow(df_))
  print("max match")
  print(sum(df_$max_match)/nrow(df_))

  # write.csv(df_, paste0(out_dir, "accuracy.csv"))

  df_min_match <- df_ %>%
    group_by(main_label, data_source) %>%
    summarise_at(vars(min_match), list(perc_min_match = mean)) %>%
    arrange(desc(perc_min_match)) 

  df_max_match <- df_ %>%
    group_by(main_label, data_source) %>%
    summarise_at(vars(max_match), list(perc_max_match = mean)) %>%
    arrange(desc(perc_max_match)) 

  #merge 
  df_matches <- df_min_match %>% 
    right_join(df_max_match, by=c("main_label","data_source")) %>%
    as.data.frame() %>%
    arrange(main_label) %>%
    mutate(label_source = paste(main_label, data_source, sep = " ")) 

  # print(df_matches)

  df_matches_ <- df_matches %>%
    select(perc_min_match, perc_max_match,label_source  ) %>%
    gather("type", "percent", -label_source )

  # print(df_matches_)

  s <- ggplot(df_matches_, aes(y =percent, x= label_source, fill=type) )+
        geom_bar(stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      labs (title = paste0("method: ",name))
  
  ggsave(plot=s, filename=paste0(out_dir, "_acc_plot.png"))
}


gsva_heatmap <- function(df_enr, out_path){
# colors from https://sashamaps.net/docs/resources/20-colors/

  # ls_col = c(
  #     "B cells memory" = "#469990", #Teal                           
  #     "B cells naive"  = "#000075", #Navy

  #     "Eosinophils" = "#42d4f4",
  #     "Neutrophils" = "#4363d8",

  #     "Monocytes" = "#800000", #Black
  #     "Macrophages M0" = "#9A6324",           
  #     "Macrophages M1"  = "#808000", #Maroon

  #     "Dendritic cells resting" = "#f58231",
  #     "Dendritic cells activated" = "#e6194B",
                          
  #     "NK cells resting" = "#ffe119",

  #     "T cells gamma delta"   = "#fabed4", #pink
  #     "T cells follicular helper" = "#ffd8b1", #apricot
  #     "T cells regulatory"  = "#fffac8", #beige
  #     "T cells CD4 naive"   = "#aaffc3", #mint
  #     "T cells CD8"         = "#dcbeff" #lavender
  # )


# Red
# #e6194B
# Green
# #3cb44b
# Yellow
# #ffe119
# Blue
# #4363d8
# Orange
# #f58231
# Purple
# #911eb4
# Cyan
# #42d4f4
# Magenta
# #f032e6
# Lime
# #bfef45
# Pink
# #fabed4
# Teal
# #469990
# Lavender
# #dcbeff
# Brown
# #9A6324
# Beige
# #fffac8
# Maroon
# #800000
# Mint
# #aaffc3
# Olive
# #808000
# Apricot
# #ffd8b1
# Navy
# #000075
# Grey
# #a9a9a9
# White
# #ffffff
# Black
# #000000
# '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 

  ls_col = c(
    "B cells memory"     = "#469990", #Teal               
    "B cells naive"       = "#000075", #Navy     
    "CD4:+ T Cell"     = "#aaffc3", #mint
    "CD8:+ T Cell"   = "#dcbeff", #lavender

    "Dendritic NT"   =  "#800000", #maroon
    "Dendritic LPS-18h"     = "#e6194B", #red
    "Dendritic R837-18h"    = "#e6194B", #red
    "Dendritic R848-18h"     = "#e6194B", #red

    "Plasmacytoid DC"  = "#fabed4", #pink
    "Myeloid DC"   = "#ffd8b1",  #Apricot
    "Myeloid DC CD123+" =  "#f032e6",   #Magenta

    "Eosinophils"   = "#000000", #black
    
    "Macrophage NT"  = "#42d4f4", #cyan         
    "Macrophage LPS-18h"     = "#bfef45", #lime
    "Macrophage Pam3CSK4-18h" = "#bfef45",#lime
    "Macrophage R837-18h"    = "#bfef45", #lime
    "Macrophage R848-18h"    = "#bfef45",#lime

    "Monocyte NT"  = "#4363d8",  #blue              
    "Monocyte LPS-18h"  = "#911eb4" ,  #purple  
    "Monocyte Pam3CSK4-18h"    = "#911eb4" ,  #purple  
    "Monocyte R837-18h"    = "#911eb4" ,  #purple     
    "Monocyte R848-18h"    = "#911eb4" ,  #purple      
    "Monocytes"               = "#911eb4" ,  #purple  
    "Non-classical monocyte"  = "#9A6324", #brown
        
    "Neutrophils" = "#a9a9a9",  #Grey
    "NK cells resting" = "#fffac8" #beige

  )

  ha <- HeatmapAnnotation(
      df = metadata, 
      name = "main_label", 
      col = list("main_label"= ls_col, "data_source" = c("Choi" = "black",
                                                       "Monaco" = "darkgrey" , 
                                                       "Song" = "lightgrey")),
       na_col = "grey",
      annotation_legend_param = list(),
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

  ###############################################################################
  # clustered unscaled heatmap 
  png(file=paste0(out_path,"_heatmap_unscaled_clustered.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht <- ComplexHeatmap::Heatmap(df_enr,
                                # name="==",
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

  draw(ht, merge_legend = TRUE)
  dev.off()


  ##############################################################################
  # clustered scaled heatmap 
  df_enr_scale <- t(scale(t(df_enr))) # scale and center rows

  png(file=paste0(out_path,"_heatmap_scaled_clustered.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht_scaled_clustered <- ComplexHeatmap::Heatmap(df_enr_scale,
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

  draw(ht_scaled_clustered, merge_legend = TRUE)
  dev.off()

  ##############################################################################
  # ordered scaled heatmap 
  ls_sample_order <- metadata %>% 
    rownames_to_column("Run") %>%
    arrange(main_label) %>%
    pull(Run)

  png(file=paste0(out_path,"_heatmap_scaled_ordered.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_scale,
                                name="Z-Score",
                                  column_order =ls_sample_order,
                                  row_order = order(rownames(df_enr_scale)), 
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

  draw(ht_order, merge_legend = TRUE)
  dev.off()

  ###############################################################################
  # ordered unscaled heatmap 
  png(file=paste0(out_path,"_heatmap_unscaled_ordered.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht <- ComplexHeatmap::Heatmap(df_enr,
                                # name="==",
                                  column_order =ls_sample_order,
                                  show_row_names= TRUE,
                                  row_order = order(rownames(df_enr_scale)), 
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

  draw(ht, merge_legend = TRUE)
  dev.off()

}


# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "GSVA"),

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
  arrange(Run) %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Read in and sort metadata 
metadata <- metadata %>% 
    arrange(Run) %>%
    tibble::column_to_rownames("Run") %>%
    select(main_label, group_label, data_source)

# Read in and sort enrichment output 
df_enr <- read.csv(file = opt$i, stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("name")
df_enr <- df_enr[ , order(names(df_enr))]

# Make heatmap of enrichment output 
gsva_heatmap(df_enr, paste0(opt$out_dir,"gsva") )

# Add row for which set has max and min score in each sample 
max <- rownames(df_enr)[apply(df_enr,2,which.max)]
df_enr["max",] <- max
min <- rownames(df_enr)[apply(df_enr,2,which.min)]
df_enr["min",] <- min

df_enr["main_label", ] <- df_sample_annotations["main_label",]
df_enr["group_label", ] <- df_sample_annotations["group_label",]
df_enr["data_source", ] <- df_sample_annotations["data_source",]

# Transpose to compare label to max and min 
# check if max value is label + "UP"
# check if min value is label + "DN"
df_ <- df_enr %>%
    t() %>%
    as.data.frame() %>%
    select(main_label,group_label, data_source, min, max) %>%
    mutate(main_label_str = gsub(" ", "_", main_label)) %>%
    mutate(group_label_str = gsub(" ", "_", group_label)) %>%
    mutate(max_main_match= ifelse((startsWith(as.character(max),
                                        as.character(main_label_str)) & 
                                endsWith(as.character(max),
                                        as.character("UP"))) , 1, 0)) %>%
    mutate(min_main_match= ifelse((startsWith(as.character(min),
                                        as.character(main_label_str)) & 
                                endsWith(as.character(min),
                                        as.character("DN"))) , 1, 0)) %>%
    mutate(max_group_match= ifelse((startsWith(as.character(max),
                                        as.character(group_label_str)) & 
                                endsWith(as.character(max),
                                        as.character("UP"))) , 1, 0)) %>%
    mutate(min_group_match= ifelse((startsWith(as.character(min),
                                        as.character(group_label_str)) & 
                                endsWith(as.character(min),
                                        as.character("DN"))) , 1, 0))  %>%
    mutate(min_main_or_group_match = ifelse((min_main_match ==1| min_group_match==1)
                                             , 1, 0))      %>%
    mutate(max_main_or_group_match = ifelse((max_main_match ==1| max_group_match==1)
                                             , 1, 0))           

cat("\n")
cat("Min exact match with main DN label: \n")
print(sum(df_$min_main_match)/nrow(df_))
cat("Max exact match with main UP label: \n")
print(sum(df_$max_main_match)/nrow(df_))

cat("\n")
cat("Min exact match with group DN label: \n")
print(sum(df_$min_group_match)/nrow(df_))
cat("Max exact match with group UP label: \n")
print(sum(df_$max_group_match)/nrow(df_))

cat("\n")
cat("Min match with main or group DN label: \n")
print(sum(df_$min_main_or_group_match)/nrow(df_))
cat("Max match with main or group UP label: \n")
print(sum(df_$max_main_or_group_match)/nrow(df_))


write.csv(df_, file = paste0(opt$o, "gsva_eval_table.csv"), row.names = FALSE)
