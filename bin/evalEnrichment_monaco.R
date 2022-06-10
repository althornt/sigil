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
    "CD4:+_T_Cell_main_DN",
    "CD4:+_T_Cell_main_UP",
    "CD4:+_T_Cell_within_DN",
    "CD4:+_T_Cell_within_UP",
    "CD8:+_T_Cell_main_DN",
    "CD8:+_T_Cell_main_UP",
    "CD8:+_T_Cell_within_DN",
    "CD8:+_T_Cell_within_UP",
    "T_cells_group_DN",
    "T_cells_group_UP"
)


DC_sets <- list(
    "Myeloid_DC_CD123+_main_DN",
    "Myeloid_DC_CD123+_main_UP",
    "Myeloid_DC_CD123+_within_DN",
    "Myeloid_DC_CD123+_within_UP",
    "Myeloid_DC_main_DN",
    "Myeloid_DC_main_UP",
    "Myeloid_DC_within_DN",
    "Myeloid_DC_within_UP",
    "Plasmacytoid_DC_main_DN",
    "Plasmacytoid_DC_main_UP",
    "Plasmacytoid_DC_within_DN",
    "Plasmacytoid_DC_within_UP",
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
    "Eosinophils_group_DN",
    "Eosinophils_group_UP",
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
monocyte_sets <- list(
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
    "Neutrophils_main_UP",
    "Neutrophils_group_DN",
    "Neutrophils_group_UP"
)

###############################
# Monoco cell types

monaco_tcells <- list(
    "Naive CD4 T cells", 
    "Naive CD8 T cells",
    "Central memory CD8 T cell",
    "Effector memory CD8 T cells",
    "Terminal effector CD8 T cells",
    "Follicular helper T cells",
    "T regulatory cells",
    "Th1 cells Th7 cells",
    "Th7 cells",
    "Th1 cells",
    "Th1/Th17 cells",
    "Th17 cells",
    "Th2 cells",
    "MAIT cells",
    "Vd2 gd T cells",
    "Non-Vd2 gd T cells",
    "Terminal effector CD4 T cells"
)

monaco_bcells <- list(
    "Naive B cells",
    "Non-switched memory B cells",
    "Exhausted B cells",
    "Switched memory B cells"
)

monaco_monocytes <- list(
    "Classical monocytes",
    "Intermediate monocytes",
    "Non classical monocytes"
)

monaco_NKcells <- list(
    "Natural killer cells"
    )

monaco_DC <- list(
    "Plasmacytoid dendritic cells",
    "Myeloid dendritic cells"
)

monaco_neutrophil <-list(
    "Low-density neutrophils"
    )


Monaco2set <- list(
  "Tcells" = list(monaco_tcells, Tcell_sets),
  "Bcells" = list(monaco_bcells, Bcell_sets),
  "Monocytes" = list(monaco_monocytes, monocyte_sets),
  "NKcells" = list(monaco_NKcells, NK_sets),
  "DC" = list(monaco_DC, DC_sets),
  "neutrophils"    =  list(monaco_neutrophil, Neutrophil_sets)
)

# print(Monaco2set)
# print("--------")
# print(Monaco2set[[1]][1])

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

# Read in and sort enrichment output 
df_enr <- read.csv(file = opt$i, stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("name")
df_enr <- df_enr[ , order(names(df_enr))]

# # Add row for which set has max and min score in each sample 
# max <- rownames(df_enr)[apply(df_enr,2,which.max)]
# df_enr["max",] <- max
# min <- rownames(df_enr)[apply(df_enr,2,which.min)]
# df_enr["min",] <- min
# df_enr["main_label", ] <- df_sample_annotations["main_label",]
# df_enr["group_label", ] <- df_sample_annotations["group_label",]
# # df_enr["data_source", ] <- df_sample_annotations["data_source",]
# write.csv(df_enr, file = paste0(opt$o, "gsva_maxmin.csv"), row.names = TRUE)

# print(unique(metadata$group_label))

ls_groups_keep <- list("T cells", "Th cells", "B cells", 
                        "Monocytes", "NK cells", "Dendritic cells",
                        "Neutrophils")     
ls_main_drop <- list("Low-density basophils","PBMCs",
                      "Plasmablasts","Progenitor cells")


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


checktop <-function(str_main_label, ls_correct_sets, df_enr){

   # Get list of samples of this type
    ls_samples <- metadata %>%
      tibble::rownames_to_column("Run")  %>%
      filter(main_label == str_main_label) %>%
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

# For each cell type use Monaco2set list to check if the max UP or min DN
# is in the correct splice set
for ( i in unique(metadata$main_label)){
  # Tcells 
  if (i %in% unlist(Monaco2set[[1]][1])){
    df_enr <- checktop(i, unlist(Monaco2set[[1]][2]), df_enr)
  } 
  # Bcells
  else if (i %in% unlist(Monaco2set[[2]][1])){
    df_enr <- checktop(i, unlist(Monaco2set[[2]][2]), df_enr)
  } 
  # Monocytes
  else if (i %in% unlist(Monaco2set[[3]][1])){
    df_enr <- checktop(i, unlist(Monaco2set[[3]][2]), df_enr)
  }
  # NKcells
  else if (i %in% unlist(Monaco2set[[4]][1])){
    df_enr <- checktop(i, unlist(Monaco2set[[4]][2]), df_enr)
  }
  # DC 
  else if (i %in% unlist(Monaco2set[[5]][1])){
    df_enr <- checktop(i, unlist(Monaco2set[[5]][2]), df_enr)
  }
  # Neutrophils
  else if (i %in% unlist(Monaco2set[[6]][1])){
    df_enr <- checktop(i, unlist(Monaco2set[[6]][2]), df_enr)
  }

  } 



df_enr["main_label", ] <- df_sample_annotations["main_label",]
df_enr["group_label", ] <- df_sample_annotations["group_label",]

write.csv(df_enr, file = paste0(opt$o, "gsva_maxmin_main.csv"), row.names = TRUE)


# Drop the samples that didnt have matching sets 
# df_enr<- df_enr %>%
#   select_if(~ !any(is.na(.)))

df_enr_main_acc <- df_enr %>%
  t() %>%
  as.data.frame()  %>%
  filter(!main_label %in% ls_main_drop) %>%
  mutate_at(c('max_up_set_or_min_down_set_match'), as.character) %>%
  mutate_at(c('max_up_set_or_min_down_set_match'), as.numeric) 

# Overall accuracy from all samples
mean(df_enr_main_acc$max_up_set_or_min_down_set_match)

# Accuracy by main label 
df_enr_main_acc_main <- df_enr_main_acc %>%
  select(max_up_set_or_min_down_set_match, main_label) %>%
  group_by(main_label ) %>%
  dplyr::summarize(Mean = mean(max_up_set_or_min_down_set_match, na.rm=FALSE)) %>%
  as.data.frame()  %>%
  arrange(desc(Mean)) 

# Accuracy by group label 
df_enr_main_acc_group <- df_enr_main_acc %>%
  select(max_up_set_or_min_down_set_match, group_label) %>%
  group_by(group_label ) %>%
  dplyr::summarize(Mean = mean(max_up_set_or_min_down_set_match, na.rm=FALSE)) %>%
  as.data.frame()  %>%
  arrange(desc(Mean)) 

print(df_enr_main_acc_main)
print(df_enr_main_acc_group)


#######################
# median heatmap by group
#########################

# group enr by the group label 
df_enr_median <- df_enr %>%
  t() %>%
  as.data.frame()  %>%
  select(-max_up_set_or_min_down_set_match,-max_up_set,-min_down_set, -main_label ) %>%
  filter(group_label %in% ls_groups_keep) %>%
  mutate_at(vars(ends_with("_UP")), funs(as.numeric(as.character(.)))) %>%
  mutate_at(vars(ends_with("_DN")), funs(as.numeric(as.character(.)))) %>%
  group_by(group_label) %>%
  summarise_all(.funs = c(median="median"))%>%
  # t() %>%
  as.data.frame()  

print(head(df_enr_median))



##############################################################################
# Heatmap all 

df_enr_median_heat <- df_enr_median %>%
  tibble::column_to_rownames("group_label") 

df_enr_median_heat[] <- sapply(df_enr_median_heat, as.numeric)
  
png(file=paste0(opt$out_dir,"/median_heatmap_group.png"),
    width = 9,
    height    = 10,
    units     = "in",
    res       = 1200)

ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat,
                              # name="Z-Score",
                                column_order =order(colnames(df_enr_median_heat)),
                                # row_order = order(rownames(df_enr_scale)), 
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


##############################################################################
# Heatmap DN

df_enr_median_heat_DN <- df_enr_median_heat %>%
    select(matches("DN_median"))


print(head(df_enr_median_heat_DN))

  png(file=paste0(opt$out_dir,"/median_heatmap_DN_group.png"),
      width = 9,
      height    = 7,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat_DN,
                                # name="Z-Score",
                                  column_order =order(colnames(df_enr_median_heat_DN)),
                                  row_order = order(rownames(df_enr_median_heat_DN)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =8),
                                  column_names_gp = grid::gpar(fontsize =7),
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

##############################################################################
# Heatmap DN

df_enr_median_heat_UP <- df_enr_median_heat %>%
    select(matches("UP_median"))


print(head(df_enr_median_heat_UP))

  png(file=paste0(opt$out_dir,"/median_heatmap_UP_group.png"),
      width = 9,
      height    = 7,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat_UP,
                                # name="Z-Score",
                                  column_order =order(colnames(df_enr_median_heat_UP)),
                                  row_order = order(rownames(df_enr_median_heat_UP)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =8),
                                  column_names_gp = grid::gpar(fontsize =7),
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

#######################
# median heatmap by main
#########################

# group enr by the group label 
df_enr_median <- df_enr %>%
  t() %>%
  as.data.frame()  %>%
  select(-max_up_set_or_min_down_set_match,-max_up_set,-min_down_set, -group_label ) %>%
  # filter(main_label %in% ls_groups_keep) %>%
  mutate_at(vars(ends_with("_UP")), funs(as.numeric(as.character(.)))) %>%
  mutate_at(vars(ends_with("_DN")), funs(as.numeric(as.character(.)))) %>%
  group_by(main_label) %>%
  summarise_all(.funs = c(median="median"))%>%
  # t() %>%
  as.data.frame()  

print(head(df_enr_median))



##############################################################################
# Heatmap all 

df_enr_median_heat <- df_enr_median %>%
  tibble::column_to_rownames("main_label") 

df_enr_median_heat[] <- sapply(df_enr_median_heat, as.numeric)
  
png(file=paste0(opt$out_dir,"/median_heatmap_main.png"),
    width = 9,
    height    = 10,
    units     = "in",
    res       = 1200)

ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat,
                              # name="Z-Score",
                                # column_order =order(colnames(df_enr_median_heat)),
                                # row_order = order(rownames(df_enr_scale)), 
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


##############################################################################
# Heatmap DN

df_enr_median_heat_DN <- df_enr_median_heat %>%
    select(matches("DN_median"))


print(head(df_enr_median_heat_DN))

  png(file=paste0(opt$out_dir,"/median_heatmap_DN_main.png"),
      width = 9,
      height    = 7,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat_DN,
                                # name="Z-Score",
                                  # column_order =order(colnames(df_enr_median_heat_DN)),
                                  # row_order = order(rownames(df_enr_median_heat_DN)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =8),
                                  column_names_gp = grid::gpar(fontsize =7),
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

##############################################################################
# Heatmap DN

df_enr_median_heat_UP <- df_enr_median_heat %>%
    select(matches("UP_median"))


print(head(df_enr_median_heat_UP))

  png(file=paste0(opt$out_dir,"/median_heatmap_UP_main.png"),
      width = 9,
      height    = 7,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat_UP,
                                # name="Z-Score",
                                  # column_order =order(colnames(df_enr_median_heat_UP)),
                                  # row_order = order(rownames(df_enr_median_heat_UP)), 
                                  show_row_names= TRUE,
                                  show_column_names = TRUE,
                                  row_names_gp = grid::gpar(fontsize =8),
                                  column_names_gp = grid::gpar(fontsize =7),
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