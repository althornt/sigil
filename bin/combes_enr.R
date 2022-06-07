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
    summarise_at(vars(min_match), list(perc_min_match = median)) %>%
    arrange(desc(perc_min_match)) 

  df_max_match <- df_ %>%
    group_by(main_label, data_source) %>%
    summarise_at(vars(max_match), list(perc_max_match = median)) %>%
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





# Arguments
option_list <- list(
  optparse::make_option(
    c("-e", "--enr"),
    type = "character",
    default = NULL,
    help = "gsva output csv"),

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
  dplyr::select(Run,main_label,group_label) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
# all_PS_meta <- rbind(all_PS, df_sample_annotations)

# Make output directories
ls_out_paths <- list("/gsva" )
for (path in ls_out_paths){

  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}
}


metadata <- metadata %>% 
    arrange(Run) %>%
    tibble::column_to_rownames("Run") %>%
    select(main_label)


# Make heatmap of ssGSEA res
# ssgsea_heatmap(dat.gct, opt$out_dir)

# Calculate and plot accuracty of ssGSEA res
# ssgsea_acc(dat.gct, opt$out_dir)


gsva_heatmap <- function(df_enr, out_path){
# colors from https://sashamaps.net/docs/resources/20-colors/

  ls_col = c(
      "Treg lung" = "#469990", #Teal                           
      "Treg kidney"  = "#000075", #Navy

      "Tcell lung" = "#42d4f4",
      "Tcell kidney" = "#4363d8"
  )

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

df_enr <- read.csv(file = opt$enr, stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("name")
df_enr <- df_enr[ , order(names(df_enr))]

# # order cols to match metadata
df_enr <- df_enr[ , order(names(df_enr))]
# gsva_heatmap(df_enr, paste0(opt$out_dir,"gsva") )

print(df_enr)
print(length(rownames(df_enr)))
print(metadata)

# Add row for which set has max and min score in each sample 
max <- rownames(df_enr)[apply(df_enr,2,which.max)]
df_enr["max",] <- max
min <- rownames(df_enr)[apply(df_enr,2,which.min)]
df_enr["min",] <- min

df_enr["main_label", ] <- df_sample_annotations["main_label",]

write.csv(df_enr, paste0(opt$out_dir,"/enr_min_max.csv"))


# Calc freq of sets with min score in T cells 
df_enr_min_tcell <- df_enr %>%
  t() %>%
  as.data.frame() %>%
  filter(main_label %in% c("Tcell kidney", "Tcell lung")) %>%
  select(main_label, min) %>%
  # count(main_label, min) %>%
  count(min) %>%
  as.data.frame() %>%
  mutate(freq = round(n / sum(n), 3)) %>% 
  arrange(desc(freq)) 
  # %>%
  # filter(grepl("DN",min))

cat("\n Tcell sets with minimum score: \n")
print(df_enr_min_tcell)

# Calc freq of sets with min score in T regs 

df_enr_min_treg <- df_enr %>%
  t() %>%
  as.data.frame() %>%
  filter(main_label %in% c("Treg kidney", "Treg lung")) %>%
  select(main_label, min) %>%
  # count(main_label, min) %>%
  count(min) %>%
  as.data.frame() %>%
  mutate(freq = round(n / sum(n), 3)) %>% 
  arrange(desc(freq)) 
  # %>%
  # filter(grepl("DN",min))
cat("\n Treg sets with minimum score: \n")
print(df_enr_min_treg)


# Calc freq of sets with max score in T cells 
df_enr_max_tcell <- df_enr %>%
  t() %>%
  as.data.frame() %>%
  filter(main_label %in% c("Tcell kidney", "Tcell lung")) %>%
  select(main_label, max) %>%
  count(max) %>%
  as.data.frame() %>%
  mutate(freq = round(n / sum(n), 3)) %>% 
  arrange(desc(freq)) 


cat("\n Tcell sets with maximum score: \n")
print(df_enr_max_tcell)

# Calc freq of sets with max score in T regs 
df_enr_max_treg <- df_enr %>%
  t() %>%
  as.data.frame() %>%
  filter(main_label %in% c("Treg kidney", "Treg lung")) %>%
  select(main_label, max) %>%
  count(max) %>%
  as.data.frame() %>%
  mutate(freq = round(n / sum(n), 3)) %>% 
  arrange(desc(freq)) 

cat("\n Treg sets with maximum score: \n")
print(df_enr_max_treg)

# print(df_enr)

# df_enr_melt <- df_enr %>%
#     select(perc_min_match, perc_max_match,label_source, -min, -max  ) %>%
#     gather("key", "val" )


# print(df_enr_melt)

print(unique(metadata$main_label))


# Calc freq of sets with max score in myeloid
df_enr_max_myeloid <- df_enr %>%
  t() %>%
  as.data.frame() %>%
  filter(main_label %in% c("myeloid lung")) %>%
  select(main_label, max) %>%
  count(max) %>%
  as.data.frame() %>%
  mutate(freq = round(n / sum(n), 3)) %>% 
  arrange(desc(freq)) 


cat("\n Myeloid sets with maximum score: \n")
print(df_enr_max_myeloid)

# Calc freq of sets with min score in myleoid
df_enr_min_myeloid <- df_enr %>%
  t() %>%
  as.data.frame() %>%
  filter(main_label %in% c("myeloid lung")) %>%
  select(main_label, min) %>%
  count(min) %>%
  as.data.frame() %>%
  mutate(freq = round(n / sum(n), 3)) %>% 
  arrange(desc(freq)) 

cat("\n Myeloid sets with min score: \n")
print(df_enr_min_myeloid)

print(head(df_enr))

##############################################################################
# Heatmap all 

  # df_enr_heat <- df_enr %>%
  #   tibble::column_to_rownames("main_label") 

df_enr[] <- sapply(df_enr, as.numeric)
print(head(df_enr))

png(file=paste0(opt$out_dir,"/gsva_heatmap.png"),
    width = 9,
    height    = 10,
    units     = "in",
    res       = 1200)

ht_order <- ComplexHeatmap::Heatmap(df_enr,
                              # name="Z-Score",
                                # column_order =order(colnames(df_enr)),
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

quit()

#######################################################

# Group input samples by main_label
# get median of each score in each group 
# rank groups by median 
df_enr_median <- df_enr %>%
  t() %>%
  as.data.frame()  %>%
  select(-min,-max ) %>%
  mutate_at(vars(ends_with("_UP")), funs(as.numeric(as.character(.)))) %>%
  mutate_at(vars(ends_with("_DN")), funs(as.numeric(as.character(.)))) %>%
  group_by(main_label) %>%
  summarise_all(.funs = c(median="median"))%>%
  # t() %>%
  as.data.frame()  




print(head(df_enr_median))
# print(df_enr_median)


# df_enr_median_melt <- df_enr_median %>%
#     # select(-main_label  ) %>%
#     melt(id = "main_label" )

# print(head(df_enr_median_melt))

# print(dim(df_enr_median_melt))
# df_enr_median_melt_UP <- df_enr_median_melt %>%
#   filter(grepl('UP_median', variable)) %>%
#   filter(variable %in% c("T_cells_group_UP_median","Monocyte_R848-18h_main_UP_median","Macrophage_LPS-18h_main_UP"))

# print(dim(df_enr_median_melt_UP))

# ggplot(data = df_enr_median_melt_UP, aes(x = main_label, y = value)) + 
#        # `geom_col()` uses `stat_identity()`: it leaves the data as is.
#        geom_col(position = 'dodge')+
#       facet_grid(~variable) +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# ggsave(filename=paste0(opt$out_dir,"/medianplot.png"))

##############################################################################
# Heatmap all 

  df_enr_median_heat <- df_enr_median %>%
    tibble::column_to_rownames("main_label") 

  df_enr_median_heat[] <- sapply(df_enr_median_heat, as.numeric)
    
  png(file=paste0(opt$out_dir,"/median_heatmap.png"),
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

  png(file=paste0(opt$out_dir,"/median_heatmap_DN.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat_DN,
                                # name="Z-Score",
                                  column_order =order(colnames(df_enr_median_heat_DN)),
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

df_enr_median_heat_UP <- df_enr_median_heat %>%
    select(matches("UP_median"))


print(head(df_enr_median_heat_UP))

  png(file=paste0(opt$out_dir,"/median_heatmap_UP.png"),
      width = 9,
      height    = 10,
      units     = "in",
      res       = 1200)

  ht_order <- ComplexHeatmap::Heatmap(df_enr_median_heat_UP,
                                # name="Z-Score",
                                  column_order =order(colnames(df_enr_median_heat_UP)),
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