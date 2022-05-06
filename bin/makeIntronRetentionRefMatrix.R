#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(uwot)
library(RColorBrewer)
library(foreach)
library(doParallel)
# library(import)

cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

# Note - the main difference between this script and the PS script is that
# For intron retention you need to consider all junctions even those that do not have 
# mutually exclusive junctions. 

make_pheatmap <- function(ls_events, label, df_meta, df_PS){
  #' Make heatmap using the pheatmap package using the given events
  #' and data
  #'
  #' @param ls_events - list of events of interest to use in heatmap
  #' @param label - label to use in output file path
  #' @param df_meta - df of metadata with Run, val, and data_source, main_label, group_label
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
      dplyr::mutate(gglabel = case_when(((p.value<.00001) | ((p.value < .01  ) & (abs(delta) > .2 )))~ get(col4label)))

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
            vjust="inward",hjust="inward",
            # nudge_x = 0.05, nudge_y = 0.05,
            check_overlap =F, col = "darkgreen", size = 2
          )

    }

    ggsave(plot = p, filename = paste0(plot_out_dir,cell_type,
                                      tag,"_",col4label,"_volcano.png"),
                                      width = 12, height = 7)

  }

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

filter_top_junction <-  function(css_df){
  # get intersection of top 1000 junctions by pvalue and the ls_clusters 
  top_css_df <- css_df  %>%
      dplyr::filter(p.value <= .01 ) %>%
      dplyr::filter(abs(delta) > .2 ) %>%
      head(2000) 

  # Note fixed this IR script to look for best junc per cluster from *all clusters*
  # But can make this faster by reducing the ls_clusters I'm looping through below to 
  # only clusters that contain a cluster from the top css_df above
  ls_keep_junctions <- list()
  for (c in ls_clusters) {
      # This base R version is much faster than using dplyr
      # Find junction with lowest p.value in the given cluster
      css_df_events <- top_css_df[top_css_df$event  %in% as.list(c),]
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
  input_label <-  substr(filename, 1, nchar(filename)-4)
  print(input_label)
  # Make output directories
  if (!dir.exists(paste0(plot_out_dir,input_label))){
    dir.create(paste0(plot_out_dir,input_label),
     recursive = TRUE, showWarnings = TRUE)
  }

  # Open mesa css
  df_css <- read.table(
            file = paste0(css_dir, filename),
            sep="\t", header = TRUE)

  # Make volcano plot before filtering
  volcano(df_css, paste0(plot_out_dir,"volcanos/"),
              input_label,"_all_junctions", list())

  # Filter to most significant junction per cluster
  df_filtered <- filter_top_junction(df_css)
  write.table(df_filtered,
              file=paste0(plot_out_dir,input_label,"_css_output_top_junc_per_cluster.tsv"),
              sep = "\t",row.names = FALSE, quote=F)

  if (nrow(df_filtered) > 0){

    # Get top negative delta events
    top_sig_by_pval_negdelta <- df_filtered %>%
        dplyr::filter(p.value <= .01 ) %>%
        dplyr::filter(delta < -.2 ) %>%
        dplyr::arrange(p.value) %>%
        head(topN) %>%
        select(event,overlapping)

    ls_top_sig_by_pval_negdelta <- top_sig_by_pval_negdelta$event

    # Make plots for top negative events
    lapply(ls_top_sig_by_pval_negdelta,  plot_event, cell_type = input_label,
          LM_type=meta_col, out_dir = plot_out_dir)

    # Get top positive delta events
    top_sig_by_pval_posdelta <- df_filtered %>%
      dplyr::filter(p.value <= .01) %>%
      dplyr::filter(delta > .2 ) %>%
      dplyr::arrange(p.value) %>%
      head(topN) %>%
      select(event,overlapping)

    ls_top_sig_by_pval_posdelta <- top_sig_by_pval_posdelta$event
    # Make plots for top positive events
    lapply(ls_top_sig_by_pval_posdelta,  plot_event, cell_type = input_label,
            LM_type=meta_col,
            out_dir = plot_out_dir )

    # Combine lists
    top_sig_neg_and_pos <- unlist(list(ls_top_sig_by_pval_negdelta,
                                      ls_top_sig_by_pval_posdelta))

    # Make volcano plot labeling top_neg_and_pos
    volcano(df_css, paste0(plot_out_dir,"volcanos/"),
            input_label,"_filtered_junctions", top_sig_neg_and_pos)
    df_top_sig_neg_and_pos <- rbind(top_sig_by_pval_negdelta, top_sig_by_pval_posdelta)

    df_top_sig_neg_and_pos$group <- comparison_label
    df_top_sig_neg_and_pos$cell_type <- input_label

    return(df_top_sig_neg_and_pos)
  }
}


import_mesa_css_within<- function(ls_cell_types, topN,  label, css_dir, meta_col){
  #' Import results from MESA compare_sample_sets runs within a broader cell type
  #' (e.g. within T-cells). Find the top signficant events, make event level
  #' plots, make heatmaps of the events. This function calls import_mesa_css(),
  #' unpack_import_css_res(), and make_pheatmap()
  #'
  #' @param ls_cell_types - list of main label cell-types that were compared in
  #' compareWithinType.R script
  #' @param topN - integer; how many of the top splicing events to use
  #' @param label - string to use to represent cell type in output files
  #' @return ls_top_events - list containing 3 list - top positive events, top
  #' negative events, and top negative and positive

  print("ls_cell_types")
  print(ls_cell_types)
  print("label")
  print(label)
  print("css_dir")
  print(css_dir)
  print("meta_col")
  print(meta_col)

  # Get output files from compareWithinType script
  ls_css_file_names <- list.files(css_dir,pattern = ".tsv")
  print("ls_css_file_names")
  print(ls_css_file_names)

  stopifnot(length(ls_css_file_names) >0)

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
  make_pheatmap(ls_top_events, paste0(label, "_diff_IR_heatmap"),
          df_metadata_subset, df_all_PS )

  # Make heatmap with this cell types events and ALL samples
  make_pheatmap(ls_top_events, paste0(label,"_diff_IR_heatmap_all_samples"),
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
metadata = read.csv(file = opt$metadata, stringsAsFactors=FALSE)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
all_PS_meta <- rbind(all_PS, df_sample_annotations)

df_clusters <- read.table(file = opt$mesa_cluster, sep="\t", header = FALSE)
df_clusters$V2 <-  strsplit(as.character(df_clusters$V2), ",")

# Make list of clusters by combining first col and the list in the second col
ls_clusters <- list()
for (row in 1:nrow(df_clusters)) {
    event_main <- df_clusters[row, "V1"]
    event_others  <- df_clusters[row, "V2"]
    row_clusters <- unlist(append(as.vector(event_others),
                                      as.character(event_main)))
    ls_clusters <- append(ls_clusters, list(sort(row_clusters) ))
}

# Reduce to unique clusters; remove A:B , B:A
ls_clusters <- unique(ls_clusters)
print("list of all clusters:")
print(length(ls_clusters))
print(tail(ls_clusters))

# Make output directories
ls_out_paths <- list("/main_label/volcanos","/group_label/volcanos","/within_group/volcanos")
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
                                  pattern = ".tsv")

print(ls_main_label_css_file_names)
# Import, find signficant events and plot each one
ls_main_label_res <- foreach(i=ls_main_label_css_file_names,
                      .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_mesa_css(
      filename = i,
      topN = 20,
      meta_col="main_label",
      plot_out_dir = paste0(opt$out_dir,"/main_label/"),
      css_dir=paste0(opt$dir_group_comp,"/compare_main_label/mesa_css_outputs/"),
      comparison_label = "main_label"
      )
  }

df_main_label_res <- dplyr::bind_rows(ls_main_label_res)
print(head(df_main_label_res))
print(dim(df_main_label_res))

# Make heatmap using top events
make_pheatmap(df_main_label_res$event, "main_label_diff_IR_heatmap", metadata, all_PS )

########################################
# Import group_label 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_group_label_css_file_names <- list.files(
                                  paste0(opt$dir_group_comp,
                                  "/compare_group_label/mesa_css_outputs/"),
                                  pattern = ".tsv")

# Import, find signficant events and plot each one
ls_group_label_res <- foreach(i=ls_group_label_css_file_names,
                      .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_mesa_css(
      filename = i,
      topN = 20,
      meta_col="group_label",
      plot_out_dir = paste0(opt$out_dir,"/group_label/"),
      css_dir=paste0(opt$dir_group_comp,"/compare_group_label/mesa_css_outputs/"),
      comparison_label = "group_label"
      )
  }

df_group_label_res <- dplyr::bind_rows(ls_group_label_res)
print(head(df_group_label_res))
print(dim(df_group_label_res))

# Make heatmap using top events
make_pheatmap(df_group_label_res$event, "group_label_diff_IR_heatmap", metadata, all_PS)

########################################
# Import within comparisons
#######################################

# Match main cell type to group labels 
ls_group_cell_types <- unlist(unique(metadata[["group_label"]]))
ls_main_cell_types <- unlist(unique(metadata[["main_label"]]))

group2main <- list()
for (i in ls_group_cell_types){

  m <- metadata %>%
      filter(group_label == i) 
  if (length(unique(m$main_label)) > 1) {
      group2main<- append(group2main,list(list(i, unique(m$main_label)) ))
  }
}

# Run in parallel
ls_within_res <- foreach(i=group2main, .packages=c('magrittr','dplyr','ggplot2','pheatmap')
  )%dopar%{
    import_mesa_css_within(
      ls_cell_types = i[[2]] , 
      label = i[[1]], 
      topN = 20,
      css_dir = paste0(opt$dir_within_comp, "/mesa_css_outputs/"),
      meta_col = "main_label")
      }

df_within_res <- dplyr::bind_rows(ls_within_res, .id = "column_label")
print(df_within_res)

# Make heatmap using top events
make_pheatmap(df_within_res$event, "withinGroup_IR_heatmap", metadata, all_PS)

########################################
# Combined reference matrix
######################################
df_combined_res <- dplyr::bind_rows(list("group_label"  = df_group_label_res,
                                        "main_label"= df_main_label_res,
                                        "within_group" = df_within_res),
                                        .id = "column_label2")
print(df_combined_res)

write.table(df_combined_res,
            file=paste0(opt$out_dir,"/combined_main_group_withingroup_combinedRefMat.tsv"),
            sep = "\t",row.names = FALSE, quote=F)

# Make heatmap using top events
make_pheatmap(df_combined_res$event, "combined_main_group_withingroup_IR_heatmap", metadata, all_PS)

