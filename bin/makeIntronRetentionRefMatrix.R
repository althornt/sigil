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

# import::from("makeSplicingRefMatrix.R")
# source("/mnt/sigil/bin/makeSplicingRefMatrix.R")
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
    dplyr::filter(event %in% ls_events) %>%
    dplyr::select(noquote(order(colnames(.)))) %>%
    tibble::column_to_rownames('event')

  for (val in list("LM22", "LM6")){
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
            nudge_x = 0.05, nudge_y = 0.05,
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
            nudge_x = 0.05, nudge_y = 0.05,
            check_overlap =F, col = "darkgreen", size = 2
          )

    }

    ggsave(plot = p, filename = paste0(plot_out_dir,cell_type,
                                      tag,"_",col4label,"_volcano.png"))

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
  LM22_type <-  substr(filename, 1, nchar(filename)-4)
  print(LM22_type)

  # Make output directories
  if (!dir.exists(paste0(plot_out_dir,LM22_type))){
    dir.create(paste0(plot_out_dir,LM22_type),
     recursive = TRUE, showWarnings = TRUE)
  }

  # Open mesa css
  df_css <- read.table(
            file = paste0(css_dir, filename),
            sep="\t", header = TRUE)

  # Make volcano plot before filtering
  volcano(df_css, paste0(plot_out_dir,"volcanos/"),
              LM22_type,"_all_junctions", list())

  # Filter to most significant junction per cluster
  # df_filtered <- filter_top_junction(df_css)
  # write.table(df_filtered,
  #             file=paste0(plot_out_dir,LM22_type,"_css_output_top_junc_per_cluster.tsv"),
  #             sep = "\t",row.names = FALSE, quote=F)

  # # Make volcano plot after filtering
  # volcano(df_filtered, paste0(plot_out_dir,"volcanos/"),
  #             LM22_type,"_filtered_junctions", list())

  # Get top negative delta events
  top_sig_by_pval_negdelta <- df_css %>%
      dplyr::filter(p.value <= .01 ) %>%
      dplyr::filter(delta < -.2 ) %>%
      dplyr::arrange(p.value) %>%
      head(topN) %>%
      select(event,overlapping)

  ls_top_sig_by_pval_negdelta <- top_sig_by_pval_negdelta$event

  # Make plots for top negative events
  lapply(ls_top_sig_by_pval_negdelta,  plot_event, cell_type = LM22_type,
        LM_type=meta_col, out_dir = plot_out_dir)

  # Get top positive delta events
  top_sig_by_pval_posdelta <- df_css %>%
    dplyr::filter(p.value <= .01) %>%
    dplyr::filter(delta > .2 ) %>%
    dplyr::arrange(p.value) %>%
    head(topN) %>%
    select(event,overlapping)

  ls_top_sig_by_pval_posdelta <- top_sig_by_pval_posdelta$event
  # Make plots for top positive events
  lapply(ls_top_sig_by_pval_posdelta,  plot_event, cell_type = LM22_type,
          LM_type=meta_col,
          out_dir = plot_out_dir )

  # Combine lists
  top_sig_neg_and_pos <- unlist(list(ls_top_sig_by_pval_negdelta,
                                    ls_top_sig_by_pval_posdelta))

  # Make volcano plot labeling top_neg_and_pos
  volcano(df_css, paste0(plot_out_dir,"volcanos/"),
          LM22_type,"_filtered_junctions", top_sig_neg_and_pos)

  df_top_sig_neg_and_pos <- rbind(top_sig_by_pval_negdelta, top_sig_by_pval_posdelta)

  df_top_sig_neg_and_pos$group <- comparison_label
  df_top_sig_neg_and_pos$cell_type <- LM22_type

  return(df_top_sig_neg_and_pos)
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
    help = ""),

   optparse::make_option(
    c("-w", "--dir_within_comp"),
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

# Open files
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
all_PS_meta <- rbind(all_PS, df_sample_annotations)

df_clusters <- read.table(file = opt$mesa_cluster, sep="\t", header = FALSE)

# Make output directories
# if (!dir.exists(paste0(opt$out_dir,"/LM22/volcanos"))){
#   dir.create(paste0(opt$out_dir,"/LM22/volcanos"),
#    recursive = TRUE, showWarnings = TRUE)
# }

# if (!dir.exists(paste0(opt$out_dir,"/LM6/volcanos"))){
#   dir.create(paste0(opt$out_dir,"/LM6/volcanos"),
#    recursive = TRUE, showWarnings = TRUE)
# }

# if (!dir.exists(paste0(opt$out_dir,"/within_type/volcanos"))){
#   dir.create(paste0(opt$out_dir,"/within_type/volcanos"),
#    recursive = TRUE, showWarnings = TRUE)
# }

# if (!dir.exists(paste0(opt$out_dir,"/UMAPs"))){
#   dir.create(paste0(opt$out_dir,"/UMAPs"),
#    recursive = TRUE, showWarnings = TRUE)
# }

ls_out_paths <- list("/LM22/volcanos","/LM6/volcanos","/within_type/volcanos", "/UMAPs" )
for (path in ls_out_paths){

  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}

}

########################################
# Import LM22 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_lm22_css_file_names <- list.files(
                                  paste0(opt$dir_group_comp,
                                  "/compare_LM22/mesa_css_outputs/"),
                                  pattern = ".tsv")
ls_lm22_css_file_names <- ls_lm22_css_file_names[!ls_lm22_css_file_names %in% c("heatmaps")]
print(ls_lm22_css_file_names)
# Import, find signficant events and plot each one
ls_lm22_res <- foreach(i=ls_lm22_css_file_names,
                      .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_mesa_css(
      filename = i,
      topN = 20,
      meta_col="LM22",
      plot_out_dir = paste0(opt$out_dir,"/LM22/"),
      css_dir=paste0(opt$dir_group_comp,"/compare_LM22/mesa_css_outputs/"),
      comparison_label = "LM22"
      )
  }

df_lm22_res <- dplyr::bind_rows(ls_lm22_res)
print(head(df_lm22_res))
print(dim(df_lm22_res))

# Make heatmap using top events
make_pheatmap(df_lm22_res$event, "LM22_diff_IR_heatmap", metadata, all_PS )