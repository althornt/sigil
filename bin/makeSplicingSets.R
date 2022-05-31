#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
# library(uwot)
# library(RColorBrewer)
# library(foreach)
library(doParallel)
# cl <- makeCluster(detectCores() - 1, outfile = "")
# registerDoParallel(cl)



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

volcano <- function(df_css, plot_out_dir, cell_type, tag, ls_point_color){

  for (col4label in list("overlapping","event")){

    # For junctions under certain thresholds label with their overlapping gene or the event
    # df_css$gg_label <- NA
    df_css <- df_css %>%
      dplyr::arrange(p.value) %>%
      dplyr::mutate(gglabel = case_when(((p.value<.00001) | ((p.value < .01  ) & (abs(delta) > .3 )))~ get(col4label)))

    # If given ls_point_color, make those points red
    df_css$point_color <- NA
    # print(ls_point_color)
    if (length(ls_point_color) > 0) {
      df_css <- df_css %>%
        dplyr::arrange(p.value) %>%
        dplyr::mutate(point_color = case_when(event %in% as.vector(ls_point_color)  ~ "Event in the reference matrix",
                                              TRUE ~ "Event not in the reference matrix"))

    # Plot with point coloring
    p <- ggplot(data=df_css, aes(x=delta, y=-log10(p.value), color = point_color )) +
          geom_point(alpha = .7  ) +
          theme_minimal() + geom_vline(xintercept=c(-0.2, 0.2), col="red") +
          geom_hline(yintercept=-log10(0.01), col="red") +
          xlim(-1.0, 1.0) +
          geom_text(
            label= df_css$gglabel,
            vjust="inward",hjust="inward",
            # nudge_x = 0.05, nudge_y = 0.05,
            check_overlap =F, col = "darkgreen", size = 2
          ) +
          scale_color_manual(name = "",
            values = c("Event in the reference matrix" = "red",
                      "Event not in the reference matrix" = "black"),
            labels = c("Event in the reference matrix",
                      "Event not in the reference matrix"))+
           theme(legend.position="bottom")

    } else {
      # Plot without point coloring
      p <- ggplot(data=df_css, aes(x=delta, y=-log10(p.value) )) +
          geom_point(alpha = .7  ) +
          theme_minimal() + geom_vline(xintercept=c(-0.1, 0.1), col="red") +
          geom_hline(yintercept=-log10(0.05), col="red") +
          xlim(-1.0, 1.0) +
          geom_text(
            label= df_css$gglabel,
            # nudge_x = 0.05, nudge_y = 0.05,
            vjust="inward",hjust="inward",
            check_overlap =F, col = "darkgreen", size = 2
          )

    }

    ggsave(plot = p, filename = paste0(plot_out_dir,cell_type,
                                      tag,"_",col4label,"_volcano.png"))

  }

}

filter_top_junction <-  function(css_df){
  
  ls_keep_junctions <- list()
  for (c in ls_clusters) {
      # This base R version is much faster than using dplyr
      # Find junction with lowest p.value in the given cluster
      css_df_events <- css_df[css_df$event  %in% as.list(c),]
      css_df_top_junc <- css_df_events[order(css_df_events$p.value),][1,]
      top_junc <- as.character(css_df_top_junc$event)
      ls_keep_junctions <- append(ls_keep_junctions, top_junc)
    }

  # Filter df to top junctions
  filt_css_df <-css_df %>%
    dplyr::filter(event %in% unique(unlist(ls_keep_junctions)))

  return(filt_css_df)
}


import_mesa_css <- function(filename, topN, plot_out_dir, css_dir, meta_col, comparison_label){
  #' Import results from MESA compare sample set script to get the top N
  #' significant events into lists
  #' @param filename -
  #' @param topN -
  #' @param plot_out_dir - path to output directory
  #' @param path tom esa compare sample set output directory
  #' @return  list of top N positive evnets, top negative events, top N
  #' negative and top N positive

  # Filename to string
  label_type <-  substr(filename, 1, nchar(filename)-4)
  print(label_type)

  # # Make output directories
  # if (!dir.exists(paste0(plot_out_dir,label_type))){
  #   dir.create(paste0(plot_out_dir,label_type),
  #    recursive = TRUE, showWarnings = TRUE)
  # }

  # Open mesa css
  df_css <- read.table(
            file = paste0(css_dir, filename),
            sep="\t", header = TRUE)

  # # Make volcano plot before filtering
  # volcano(df_css, paste0(plot_out_dir,"volcanos/"),
  #             label_type,"_all_junctions", list())

  # # Filter to most significant junction per cluster
  # df_filtered <- filter_top_junction(df_css)
  # write.table(df_filtered,
  #             file=paste0(plot_out_dir,label_type,"_css_output_top_junc_per_cluster.tsv"),
  #             sep = "\t",row.names = FALSE, quote=F)

  # # Make volcano plot after filtering
  # volcano(df_filtered, paste0(plot_out_dir,"volcanos/"),
  #             label_type,"_filtered_junctions", list())

  print("all")
  print(df_css %>% dplyr::arrange(p.value) %>% head(20))

  # Get top negative delta events
  top_sig_by_pval_negdelta <- df_css %>%
      dplyr::filter(p.value <= .01 ) %>%
      dplyr::filter(delta < -.2 ) %>%
      dplyr::arrange(p.value)  %>%
      head(topN) %>%
      select(event,overlapping)

  print(head(top_sig_by_pval_negdelta))
  print(dim(top_sig_by_pval_negdelta))
  print(top_sig_by_pval_negdelta %>% count(overlapping) %>% arrange(desc(n)))
  print("---------")

  # ls_top_sig_by_pval_negdelta <- top_sig_by_pval_negdelta$event

  # # Make plots for top negative events
  # lapply(ls_top_sig_by_pval_negdelta,  plot_event, cell_type = label_type,
  #       LM_type=meta_col, out_dir = plot_out_dir)

  # Get top positive delta events
  top_sig_by_pval_posdelta <- df_css %>%
    dplyr::filter(p.value <= .01) %>%
    dplyr::filter(delta > .2 ) %>%
    dplyr::arrange(p.value)  %>%
    head(topN) %>%
    select(event,overlapping)

  print(head(top_sig_by_pval_posdelta))
  print(dim(top_sig_by_pval_posdelta))

  print(top_sig_by_pval_posdelta %>% count(overlapping) %>% arrange(desc(n)))



  return(list(label_type, droplevels(top_sig_by_pval_posdelta$event), droplevels(top_sig_by_pval_negdelta$event )))
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
    dplyr::filter(event %in% list(paste0(sig_event), paste0(LM_type), "data_source"))%>%
    t() %>%
    as.data.frame()

  df_ <- df # copy df
  colnames(df_) <- c( "PS", paste0(LM_type), "data_source") #add column names from first row

  df_ <- df_[-1,] %>% # drop first row
          dplyr::filter(PS != "NaN") %>% #drop samples with Nan
          dplyr::filter(get(LM_type) != "") #drop samaples with out the cell type label

  df_$PS <- as.numeric(as.character(df_$PS))

  p <- ggplot( df_, aes(x = get(LM_type), y = PS, color=data_source)) +
      # geom_violin() +
      geom_jitter(position=position_jitter(0.15), alpha = 0.5, size = 2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(legend.position = "right") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "cell type")

  ggsave(plot = p, filename = paste0(out_dir,cell_type,"/",sig_event,".png"))

}

make_pheatmap <- function(ls_events, label, df_meta, df_PS){
  #' Make heatmap using the pheatmap package using the given events
  #' and data
  #'
  #' @param ls_events - list of events of interest to use in heatmap
  #' @param label - label to use in output file path
  #' @param df_meta - df of metadata with Run, val, and data_source, main_label, group_label columns
  #' @param df_PS - df of MESA all PS file

  # Filter MESA all PS file to events of interest
  df_all_PS_sig_events <- df_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_events) %>%
    dplyr::select(noquote(order(colnames(.)))) %>%
    tibble::column_to_rownames('event')

  for (val in list("main_label", "group_label")){
      if (nrow(df_all_PS_sig_events) < 50){
        rowname_on_off = "T"
      } else { rowname_on_off = "F"}

      # DF to label samples(columns) with labels
      df_sample_annotations <- df_meta %>%
        dplyr::filter(paste0(val) != "") %>%
        dplyr::select(Run, val, data_source) %>%
        dplyr::arrange(Run) %>%
        tibble::column_to_rownames("Run")

      stopifnot(rownames(df_sample_annotations) == colnames(df_all_PS_sig_events))

      heatmap_res <- pheatmap(
        main = paste0(" "),
        df_all_PS_sig_events,
        # scale = "row",
        show_rownames=get(rowname_on_off),
        show_colnames=F,
        na_col = "grey",
        annotation_col = df_sample_annotations)

      save_pheatmap_pdf(
        heatmap_res,
        paste0(opt$out_dir,"/",label,"_",val,"_rowname",rowname_on_off,".pdf"))
      }
}

import_mesa_css_within<- function(ls_cell_types, topN,  label, css_dir, meta_col){
  #' Import results from MESA compare_sample_sets runs within a broader cell type
  #' (e.g. within T-cells). Find the top signficant events, make event level
  #' plots, make heatmaps of the events. This function calls import_mesa_css(),
  #' unpack_import_css_res(), and make_pheatmap()
  #'
  #' @param ls_cell_types - list of main_label cell-types that were compared in
  #' compareWithinType.R script
  #' @param topN - integer; how many of the top splicing events to use
  #' @param label - string to use to represent cell type in output files
  #' @return ls_top_events - list containing 3 list - top positive events, top
  #' negative events, and top negative and positive


  # Get output files from compareWithinType script
  ls_css_file_names <- list.files(css_dir,pattern = ".tsv")

  # For input cell type list , convert to filename
  ls_cell_types_file <- c()
  for (val in ls_cell_types){
    new_val <- paste0(gsub(" ","_", val), ".tsv")
    ls_cell_types_file <- append(ls_cell_types_file, new_val)
  }

  # Intersect with the files that exist (Not all will have a mesa css output )
  ls_css_file_names_cell_type  <- intersect(ls_css_file_names, ls_cell_types_file)

  #Import files, find top signficant events, plot each event
  ls_res <- lapply(
                    ls_css_file_names_cell_type,
                    topN=topN,
                    import_mesa_css,
                    meta_col =meta_col,
                    plot_out_dir =  paste0(opt$out_dir,"/within_group/"),
                    css_dir =  css_dir,
                    comparison_label = "withinType")

  # If only 2 CSS files, they should have identical resuls
  # (comparing A vs B then B) , so only return the results of one
  if (length(ls_css_file_names_cell_type) <= 2) {
      #only keep first result
      df_res <- ls_res[[1]]
  } else {
      df_res <- dplyr::bind_rows(ls_res)
  }

  ls_top_events <- df_res$event

  # Filter metadata
  df_metadata_subset <- metadata %>%
    dplyr::filter(main_label %in% ls_cell_types) %>%
    droplevels(.) %>%
    dplyr::arrange(Run)

  # Filter all PS
  df_all_PS <- all_PS %>%
    dplyr::select(as.vector(unlist(df_metadata_subset$Run)))

  # Make heatmap with this cell types events only within this cell types samples
  make_pheatmap(ls_top_events, paste0(label, "_diff_splicing_heatmap"),
          df_metadata_subset, df_all_PS )

  # Make heatmap with this cell types events and ALL samples
  make_pheatmap(ls_top_events, paste0(label,"_diff_splicing_heatmap_all_samples"),
          metadata, all_PS )

  return(df_res)
}

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--mesa_PS"),
    type = "character",
    default = NULL,
    help = " `mesa_allPS.tsv` input file "),

  optparse::make_option(
    c("-c", "--mesa_cluster"),
    type = "character",
    default = NULL,
    help = "full path to mesa cluster file"),

  optparse::make_option(
    c("-g", "--dir_group_comp"),
    type = "character",
    default = NULL,
    help = "path to css outptus from runMESAcompare"),

   optparse::make_option(
    c("-w", "--dir_within_comp"),
    type = "character",
    default = NULL,
    help = "path to css outputs from compare within type"),

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

all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
all_PS_meta <- rbind(all_PS, df_sample_annotations)

# df_clusters <- read.table(file = opt$mesa_cluster, sep="\t", header = FALSE)
# # Remove rows where second column is empty (no ME Junction) and format
# df_clusters_filter <- df_clusters %>% dplyr::filter(V2 != "")
# df_clusters_filter$V2 <-  strsplit(as.character(df_clusters_filter$V2), ",")

# # Make list of clusters by combining first col and the list in the second col
# ls_clusters <- list()
# for (row in 1:nrow(df_clusters_filter)) {
#     event_main <- df_clusters_filter[row, "V1"]
#     event_others  <- df_clusters_filter[row, "V2"]
#     row_clusters <- unlist(append(as.vector(event_others),
#                                       as.character(event_main)))
#     ls_clusters <- append(ls_clusters, list(sort(row_clusters) ))
# }

# # Reduce to unique clusters; remove A:B , B:A
# ls_clusters <- unique(ls_clusters)
# print("list of all clusters:")
# print(length(ls_clusters))
# print(tail(ls_clusters))

# Make output directories
ls_out_paths <- list("/gmt" )
for (path in ls_out_paths){

  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}
}

########################################
# Import main_label 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_main_label_css_file_names <- list.files(
                                  paste0(opt$dir_group_comp,
                                  "/compare_main_label/mesa_css_outputs/"),
                                  pattern = ".tsv")[1:3]

# Import, find signficant events and plot each one
ls_main_label_res <- foreach(i=ls_main_label_css_file_names,
                      .packages=c('magrittr','dplyr','ggplot2')) %do% {
    import_mesa_css(
      filename = i,
      topN = 200,
      meta_col="main_label",
      plot_out_dir = paste0(opt$out_dir,"/main_label/"),
      css_dir=paste0(opt$dir_group_comp,"/compare_main_label/mesa_css_outputs/"),
      comparison_label = "main_label"
      )
  }

print(ls_main_label_res)

#Check its existence
if (file.exists(paste0(opt$out_dir,"/main_set.gmt"))) {
  #Delete file if it exists
  file.remove(paste0(opt$out_dir,"/main_set.gmt"))
}

outfile = file(paste0(opt$out_dir,"/gmt/main_set.gmt"), open = 'a') # open in “a”ppend mode

for (i in ls_main_label_res){
  name = i[1][1]
  print(name)
  # print(name, "\t", "desc" )

  ls_UP <- i[[2]]
  cat(paste0(name,"_UP"),"desc" ,paste0(ls_UP, sep = "\t"), '\n', file = outfile, sep = ' ', append = TRUE)

  ls_DN <- i[[3]]
  cat(paste0(name,"_DN"),"desc" ,paste0(ls_DN, sep = "\t"), '\n', file = outfile, sep = ' ', append = TRUE)

}


# ########################################
# # Import group_label 1 vs all comparisons
# #######################################

# # Get all outputs from compare sample sets 1 vs all comparisons
# ls_group_label_css_file_names <- list.files(
#                                   paste0(opt$dir_group_comp,
#                                   "/compare_group_label/mesa_css_outputs/"),
#                                   pattern = ".tsv")
# ls_group_label_css_file_names <- ls_group_label_css_file_names[!ls_group_label_css_file_names %in% c("heatmaps")]

# # Import, find signficant events and plot each one
# ls_group_label_res <- foreach(i=ls_group_label_css_file_names,
#                       .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
#     import_mesa_css(
#       filename = i,
#       topN = 20,
#       meta_col="group_label",
#       plot_out_dir = paste0(opt$out_dir,"/group_label/"),
#       css_dir=paste0(opt$dir_group_comp,"/compare_group_label/mesa_css_outputs/"),
#       comparison_label = "group_label"
#       )
#   }

# df_group_label_res <- dplyr::bind_rows(ls_group_label_res)
# print(head(df_group_label_res))
# print(dim(df_group_label_res))

# # Make heatmap using top events
# make_pheatmap(df_group_label_res$event, "group_label_diff_splicing_heatmap", metadata, all_PS)

# ########################################
# # Import within group comparisons
# #######################################
# # Match main cell type to group labels 
# ls_group_cell_types <- unlist(unique(metadata[["group_label"]]))
# ls_main_cell_types <- unlist(unique(metadata[["main_label"]]))

# group2main <- list()
# for (i in ls_group_cell_types){
#   print(i)
#   print("==")
#   m <- metadata %>%
#       filter(group_label == i) 
#   if (length(unique(m$main_label)) > 1) {
#       group2main<- append(group2main,list(list(i, unique(m$main_label)) ))
#   }
# }

# # Run in parallel
# ls_within_res <- foreach(i=group2main, .packages=c('magrittr','dplyr','ggplot2','pheatmap')
#   )%dopar%{
#     import_mesa_css_within(
#       ls_cell_types = i[[2]] , 
#       label = i[[1]], 
#       topN = 20,
#       css_dir = paste0(opt$dir_within_comp, "/mesa_css_outputs/"),
#       meta_col = "main_label")
#       }

# df_within_res <- dplyr::bind_rows(ls_within_res, .id = "column_label")
# print(df_within_res)

# # Make heatmap using top events
# make_pheatmap(df_within_res$event, "withinGroup_splicing_heatmap", metadata, all_PS)

# ########################################
# # Combined reference matrix
# ######################################
# df_combined_res <- dplyr::bind_rows(list("group_label"  = df_group_label_res,
#                                         "main_label"= df_main_label_res,
#                                         "within_group" = df_within_res),
#                                         .id = "column_label2")
# print(df_combined_res)

# write.table(df_combined_res,
#             file=paste0(opt$out_dir,"/combined_main_group_withingroup_combinedRefMat.tsv"),
#             sep = "\t",row.names = FALSE, quote=F)

# # Make heatmap using top events
# make_pheatmap(df_combined_res$event, "combined_main_group_withingroup_splicing_heatmap", metadata, all_PS)