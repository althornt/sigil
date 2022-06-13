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



Bcell_sets <- list(
    "B_cells_group_DN",
    "B_cells_group_UP",
    "B_cells_memory_main_DN",
    "B_cells_memory_main_UP",
    "B_cells_memory_within_DN",
    "B_cells_memory_within_UP",
    "B_cells_naive_main_DN",
    "B_cells_naive_main_UP",
    "B_cells_naive_within_DN",
    "B_cells_naive_within_UP"
)

Tcell_sets <- list(
   "T_Cells_CD4:+_main_DN",
   "T_Cells_CD4:+_main_UP",
   "T_Cells_CD4:+_within_DN",
   "T_Cells_CD4:+_within_UP",
   "T_Cells_CD8:+_main_DN",
   "T_Cells_CD8:+_main_UP",
   "T_Cells_CD8:+_within_DN",
   "T_Cells_CD8:+_within_UP",
   "T_cells_group_DN",
   "T_cells_group_UP"
)


DC_sets <- list(
 "DC_Myeloid_CD123+_main_DN",
 "DC_Myeloid_CD123+_main_UP",
 "DC_Myeloid_CD123+_within_DN",
 "DC_Myeloid_CD123+_within_UP",
 "DC_Myeloid_main_DN",
 "DC_Myeloid_main_UP",
 "DC_Myeloid_within_DN",
 "DC_Myeloid_within_UP",
 "DC_Plasmacytoid_main_DN",
 "DC_Plasmacytoid_main_UP",
 "DC_Plasmacytoid_within_DN",
 "DC_Plasmacytoid_within_UP",
 "Dendritic_LPS-18h_main_DN",
 "Dendritic_LPS-18h_main_UP",
 "Dendritic_LPS-18h_within_DN",
 "Dendritic_LPS-18h_within_UP",
 "Dendritic_NT_main_DN",
 "Dendritic_NT_main_UP",
 "Dendritic_NT_within_DN",
 "Dendritic_NT_within_UP",
 "Dendritic_R837-18h_main_DN",
 "Dendritic_R837-18h_main_UP",
 "Dendritic_R837-18h_within_DN",
 "Dendritic_R837-18h_within_UP",
 "Dendritic_R848-18h_main_DN",
 "Dendritic_R848-18h_main_UP",
 "Dendritic_R848-18h_within_DN",
 "Dendritic_R848-18h_within_UP",
 "Dendritic_cells_group_DN",
 "Dendritic_cells_group_UP"
)

Eos_sets <- list(
    "Eosinophils_main_DN",
    "Eosinophils_main_UP"
)

Macrophage_sets <- list(
    "Macrophage_LPS-18h_main_DN",
    "Macrophage_LPS-18h_main_UP",
    "Macrophage_LPS-18h_within_DN",
    "Macrophage_LPS-18h_within_UP",
    "Macrophage_NT_main_DN",
    "Macrophage_NT_main_UP",
    "Macrophage_NT_within_DN",
    "Macrophage_NT_within_UP",
    "Macrophage_Pam3CSK4-18h_main_DN",
    "Macrophage_Pam3CSK4-18h_main_UP",
    "Macrophage_Pam3CSK4-18h_within_DN",
    "Macrophage_Pam3CSK4-18h_within_UP",
    "Macrophage_R837-18h_main_DN",
    "Macrophage_R837-18h_main_UP",
    "Macrophage_R837-18h_within_DN",
    "Macrophage_R837-18h_within_UP",
    "Macrophage_R848-18h_main_DN",
    "Macrophage_R848-18h_main_UP",
    "Macrophage_R848-18h_within_DN",
    "Macrophage_R848-18h_within_UP",
    "Macrophages_group_DN",
    "Macrophages_group_UP"
)
Monocyte_sets <- list(
    "Monocyte_LPS-18h_main_DN",
    "Monocyte_LPS-18h_main_UP",
    "Monocyte_LPS-18h_within_DN",
    "Monocyte_LPS-18h_within_UP",
    "Monocyte_NT_main_DN",
    "Monocyte_NT_main_UP",
    "Monocyte_NT_within_DN",
    "Monocyte_NT_within_UP",
    "Monocyte_Pam3CSK4-18h_main_DN",
    "Monocyte_Pam3CSK4-18h_main_UP",
    "Monocyte_Pam3CSK4-18h_within_DN",
    "Monocyte_Pam3CSK4-18h_within_UP",
    "Monocyte_R837-18h_main_DN",
    "Monocyte_R837-18h_main_UP",
    "Monocyte_R837-18h_within_DN",
    "Monocyte_R837-18h_within_UP",
    "Monocyte_R848-18h_main_DN",
    "Monocyte_R848-18h_main_UP",
    "Monocyte_R848-18h_within_DN",
    "Monocyte_R848-18h_within_UP",
    "Monocytes_main_DN",
    "Monocytes_main_UP",
    "Monocytes_within_DN",
    "Monocytes_within_UP",
    "Non-classical_monocyte_main_DN",
    "Non-classical_monocyte_main_UP",
    "Non-classical_monocyte_within_DN",
    "Non-classical_monocyte_within_UP",
    "Monocytes_group_DN",
    "Monocytes_group_UP"
)

NK_sets <- list(
    "NK_cells_group_DN",
    "NK_cells_group_UP",
    "NK_cells_resting_main_DN",
    "NK_cells_resting_main_UP"
)

Neutrophil_sets <- list(
    "Neutrophils_main_DN",
    "Neutrophils_main_UP"
)

Combes2set <- list(
  "Tcell" = list(list("Tcell"), Tcell_sets),
  "Treg" = list(list("Treg"), Bcell_sets),
  "myeloid" = list(list("myeloid"), unlist(Monocyte_sets, Macrophage_sets,DC_sets ))
)

checktop <-function(str_main_label, ls_correct_sets, df_enr){

   # Get list of samples of this type
    ls_samples <- metadata %>%
      tibble::rownames_to_column("Run")  %>%
      filter(group_label == str_main_label) %>%
      pull(Run)
    
    # Add 1 if either max or min is a match ; 0 if neither are a match 
    for (sample in ls_samples){

      if ((df_enr["max_up_set",sample] %in% ls_correct_sets) ||
          (df_enr["min_down_set",sample] %in% ls_correct_sets)){

          df_enr["max_up_set_or_min_down_set_match", sample] <- 1

      } else {
         df_enr["max_up_set_or_min_down_set_match", sample] <- 0
      }      
    }
  return(df_enr)
}


# Arguments
option_list <- list(
  optparse::make_option(
    c("-s", "--splice_set"),
    type = "character",
    default = NULL,
    help = ""),

  optparse::make_option(
    c("-g", "--gene_set"),
    type = "character",
    default = NULL,
    help = ""),

  optparse::make_option(
    c("-i", "--immuneSigDB"),
    type = "character",
    default = NULL,
    help = ""),

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
  dplyr::select(Run,main_label,group_label) %>%
  # dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Read in and sort metadata 
metadata <- metadata %>% 
    arrange(Run) %>%
    tibble::column_to_rownames("Run") %>%
    select(main_label, group_label)
    # select(main_label, group_label, data_source)

print(unique(metadata$main_label))
print(unique(metadata$group_label))

# Read in and sort enrichment outputs 
df_gene_enr <- read.csv(file = opt$gene_set, stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("name")
df_gene_enr <- df_gene_enr[ , order(names(df_gene_enr))]

df_splice_enr <- read.csv(file = opt$splice_set, stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("name")
df_splice_enr <- df_splice_enr[ , order(names(df_splice_enr))]

df_imsig_enr <- read.csv(file = opt$immuneSigDB, stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("name")
df_imsig_enr <- df_imsig_enr[ , order(names(df_imsig_enr))]

print(dim(df_gene_enr))
print(dim(df_splice_enr))
print(dim(df_imsig_enr))

print(dim(df_imsig_enr))

calc <- function(df_enr, name){

  df_enr_UP <- df_enr %>%
      tibble::rownames_to_column("Sets") %>%
      dplyr::filter(grepl("_UP", Sets)) %>%
      tibble::column_to_rownames("Sets") 

  df_enr_DN <- df_enr %>%
      tibble::rownames_to_column("Sets") %>%
      dplyr::filter(grepl("_DN", Sets)) %>%
      tibble::column_to_rownames("Sets") 

  max <- rownames(df_enr_UP)[apply(df_enr_UP,2,which.max)]
  df_enr["max_up_set",] <- max
  min <- rownames(df_enr_DN)[apply(df_enr_DN,2,which.min)]
  df_enr["min_down_set",] <- min

  df_enr["group_label", ] <- df_sample_annotations["group_label",]

  # For each cell type use list to check if the max UP or min DN
  # is in the correct splice set
  for ( i in unique(metadata$group_label)){

    # Tcells 
    if (i %in% unlist(Combes2set[[1]][1])){
      df_enr <- checktop(i, unlist(Combes2set[[1]][2]), df_enr)
    } 
    # Treg
    else if (i %in% unlist(Combes2set[[2]][1])){
      df_enr <- checktop(i, unlist(Combes2set[[2]][2]), df_enr)
    } 
    # Myeloid
    else if (i %in% unlist(Combes2set[[3]][1])){
      df_enr <- checktop(i, unlist(Combes2set[[3]][2]), df_enr)
    }
    } 

    df_enr_main_acc <- df_enr %>%
      t() %>%
      as.data.frame()  %>%
      filter(group_label %in% list("Tcell", "Treg", "myeloid")) %>%
      mutate_at(c('max_up_set_or_min_down_set_match'), as.character) %>%
      mutate_at(c('max_up_set_or_min_down_set_match'), as.numeric) 

    # Overall accuracy from all samples
    overall_mean <- mean(df_enr_main_acc$max_up_set_or_min_down_set_match)

    # Accuracy by group label 
    df_enr_main_acc_group <- df_enr_main_acc %>%
      select(max_up_set_or_min_down_set_match, group_label) %>%
      group_by(group_label ) %>%
      dplyr::summarize(Mean = mean(max_up_set_or_min_down_set_match, na.rm=FALSE)) %>%
      as.data.frame()  %>%
      arrange(desc(Mean)) 

  print(df_enr_main_acc_group)

  #################
  # heatmap all samples 

  df_enr_heat <-  df_enr %>%
      t() %>%
      as.data.frame()  %>%
      select(-max_up_set_or_min_down_set_match,-max_up_set,-min_down_set ) %>%
      arrange(rownames(.)) %>%
      filter(group_label %in% list("Tcell", "Treg", "myeloid"))  %>%
      select(-group_label) %>%
      mutate_at(vars(ends_with("_UP")), funs(as.numeric(as.character(.)))) %>%
      mutate_at(vars(ends_with("_DN")), funs(as.numeric(as.character(.)))) 

  png(file=paste0(opt$out_dir,"/", name, "_heatmap_group.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  metadata_harow <- metadata %>% 
        select(group_label) %>%
        filter(group_label %in% list("Tcell", "Treg", "myeloid"))

  ha_row <- rowAnnotation(
      df = metadata_harow, 
      name = "label", 
      col = list( "group_label" = c("myeloid" = "black",
                                  "Tcell" = "skyblue" , 
                                  "Treg" = "lightgrey")
                                  # "live" = "pink", 
                                  # "epcam" = "green",
                                  # "tumor" = "purple")
                                  ),
      show_legend = TRUE,
      gp = gpar(col = NA)
      )
  ht_order <- ComplexHeatmap::Heatmap(df_enr_heat,
                                # name="Z-Score",
                                  column_order =order(colnames(df_enr_heat)),
                                  # row_order = order(rownames(df_enr_heat)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =5),
                                  column_names_gp = grid::gpar(fontsize =5),
                                  right_annotation=ha_row,
                                  # heatmap_legend_param = list(legend_gp = gpar(fontsize = 3)),
                                  # show_row_dend = TRUE,
                                  heatmap_legend_param = list(
                                  # legend_direction = "horizontal", 
                                  # legend_height = unit(1, "cm"),
                                  legend_gp = gpar(fontsize = 5))
                )

  draw(ht_order, merge_legend = TRUE)
  dev.off()

  ###############
  # median
  #################

  # group enr by the group label 
  df_enr_median <- df_enr %>%
    t() %>%
    as.data.frame()  %>%
    select(-max_up_set_or_min_down_set_match,-max_up_set,-min_down_set ) %>%
    filter(group_label %in% list("Tcell", "Treg", "myeloid"))  %>%
    mutate_at(vars(ends_with("_UP")), funs(as.numeric(as.character(.)))) %>%
    mutate_at(vars(ends_with("_DN")), funs(as.numeric(as.character(.)))) %>%
    group_by(group_label) %>%
    summarise_all(.funs = c(median="median"))%>%
    # t() %>%
    as.data.frame()  

  df_enr_median_heat <- df_enr_median %>%
    tibble::column_to_rownames("group_label") 

  df_enr_median_heat[] <- sapply(df_enr_median_heat, as.numeric)

  # Make DN df
  df_enr_median_heat_DN <- df_enr_median_heat %>%
      select(matches("DN_median"))
  names(df_enr_median_heat_DN) = gsub(pattern = "_DN_median", replacement = "",
                                      x = names(df_enr_median_heat_DN))
  # Make UP df
  df_enr_median_heat_UP <- df_enr_median_heat %>%
      select(matches("UP_median"))
  names(df_enr_median_heat_UP) = gsub(pattern = "_UP_median", replacement = "",
                                      x = names(df_enr_median_heat_UP))

  # Find sets in UP not in DN and add to DN
  DN_missing <- names(df_enr_median_heat_UP)[!(names(df_enr_median_heat_UP) %in% names(df_enr_median_heat_DN))]
  df_enr_median_heat_DN[,DN_missing] <- NA

  # Find sets in DN not in UP and add to UP
  UP_missing <- names(df_enr_median_heat_DN)[!(names(df_enr_median_heat_DN) %in% names(df_enr_median_heat_UP))]
  df_enr_median_heat_UP[,UP_missing] <- NA

  # Add missing columns with NAs
  df_enr_median_heat_UP <-df_enr_median_heat_UP[ , order(names(df_enr_median_heat_UP))]
  df_enr_median_heat_DN <- df_enr_median_heat_DN[ , order(names(df_enr_median_heat_DN))]


  ht_med_dn <- ComplexHeatmap::Heatmap(df_enr_median_heat_DN,
                                  name = "",
                                  row_title = "DN",
                                  column_order =order(colnames(df_enr_median_heat_DN)),
                                  row_order = order(rownames(df_enr_median_heat_DN)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =8),
                                  column_names_gp = grid::gpar(fontsize =7)
                                  # top_annotation=ha,
                                  # heatmap_legend_param = list(legend_gp = gpar(fontsize = 3)),
                                  # show_row_dend = TRUE,
                                  # heatmap_legend_param = list(
                                  # legend_direction = "horizontal", 
                                  # legend_height = unit(1, "cm"),
                                  # legend_gp = gpar(fontsize = 5))
                )

  ht_med_up <- ComplexHeatmap::Heatmap(df_enr_median_heat_UP,
                                  row_title = "UP",
                                  column_order =order(colnames(df_enr_median_heat_UP)),
                                  row_order = order(rownames(df_enr_median_heat_UP)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =8),
                                  column_names_gp = grid::gpar(fontsize =7),
                                  show_heatmap_legend = FALSE

                                  # top_annotation=ha,
                                  # show_row_dend = TRUE,
                                  # heatmap_legend_param = list(
                                  # legend_direction = "horizontal", 
                                  # legend_height = unit(1, "cm"),
                                  # legend_gp = gpar(fontsize = 5)
                                  )

  # combined gene and splice heatmaps
  png(file=paste0(opt$out_dir,name,"_median_heatmap_UP_DN_group.png"),
      width = 10,
      height    = 4,
      units     = "in",
      res       = 1200)

  ht_list = ht_med_up %v% ht_med_dn
  draw(ht_list, 
      merge_legend = TRUE,
      column_title = name, column_title_gp = gpar(fontsize = 16))
  dev.off()

  return(list(overall_mean,df_enr_main_acc, df_enr_main_acc_group ))

}

# Run gene
ls_gene_res <- calc(df_gene_enr, "gene")
cat("\n")
cat("\n Average gene set sensitivity:")
cat("\n", paste0(ls_gene_res[1]), "\n")
cat("\n")

# write.csv(ls_gene_res[2], file = paste0(opt$o, "gene_set_eval_by_group_label.csv"), row.names = TRUE)
# # write.csv(ls_gene_res[3], file = paste0(opt$o, "gene_set_eval_by_group_label_mean.csv"), row.names = TRUE)

# Run splice
ls_splice_res <- calc(df_splice_enr, "splice")
cat("\n")
cat("\n Average splice set sensitivity:")
cat("\n", paste0(ls_splice_res[1]), "\n")
cat("\n")

write.csv(ls_splice_res[2], file = paste0(opt$o, "splice_set_eval_by_group_label.csv"), row.names = TRUE)

# write.csv(df_enr_main_acc_main, file = paste0(opt$o, "eval_by_main_label.csv"), row.names = TRUE)
# write.csv(df_enr_main_acc_group, file = paste0(opt$o, "eval_by_group_label.csv"), row.names = TRUE)

# Merge gene and splice table
cat("\n")
df_merge <- ls_splice_res[3][[1]] %>%
    rename(Splice = Mean) %>%
    left_join(ls_gene_res[3][[1]], by = "group_label") %>%
    rename(Gene = Mean)
df_merge
cat("\n")





