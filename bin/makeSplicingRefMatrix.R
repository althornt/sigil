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
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)


# This script ...
# (1) Reads and plots the results from runMESAcompare.R , which does 1 vs all
#     comparisons using uses LM22 and LM6 cell sub types
# (2) Reads and plots the results from compareWithinType.R, which does 1 vs all
#     comparisons within T-cells, Monocytes/macrophages, etc
# (3) Makes vaarious reference matrices using events identified in the different
#     comparisons

##########################
# Functions
##########################

make_umap <- function(num_neighbor,meta_col,df_PCA,out_path) {

  set.seed(123)

  # Make color palette
  n <- length(unique(metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(df_PCA, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(df_PCA)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, metadata)

  # Plot UMAP
  gg <- ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP MESA: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste(out_path,meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)
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

volcano <- function(df_css, plot_out_dir, cell_type, tag){

  df_css$gg_label <- NA
  df_css <- df_css %>%
    dplyr::arrange() %>%
    dplyr::mutate(gglabel = case_when(((p.value<.00001) | ((p.value < .01  ) & (abs(delta) > .3 )))~ overlapping))

  p <- ggplot(data=df_css, aes(x=delta, y=-log10(p.value) )) +
        geom_point(alpha = .7) +
        theme_minimal() + geom_vline(xintercept=c(-0.1, 0.1), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red") +
        xlim(-1.0, 1.0) +
        geom_text(
          label= df_css$gglabel,
          nudge_x = 0.05, nudge_y = 0.05,
          check_overlap =F, col = "darkgreen", size = 2
        )

  ggsave(plot = p, filename = paste0(plot_out_dir,cell_type,tag,"_volcano.png"))

}

filter_top_junction <-  function(css_df){

  print("filter function.................")

  ls_keep_junctions <- list()

  for (c in ls_clusters) {
      # print("-------------------------------")
      # print(c)
      #
      #
      # print(css_df %>%
      #     dplyr::filter(event %in% c) %>%
      #     dplyr:::arrange(p.value))

      top_junc  <- css_df %>%
          dplyr::filter(event %in% c) %>%
          dplyr:::arrange(p.value) %>%
          head(1) %>%
          pull(event) %>%
          droplevels()

      # print("top_junc:")
      # print(top_junc)
      ls_keep_junctions <- append(ls_keep_junctions, as.character(top_junc))
    }

  print("final_list")
  print(length(unique(unlist(ls_keep_junctions))))
  # Filter df to top junctions
  filt_css_df <-css_df %>%
    dplyr::filter(event %in% unique(unlist(ls_keep_junctions)))
  return(filt_css_df)

}

import_mesa_css <- function(filename, topN, plot_out_dir, css_dir, meta_col){
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
  df <- read.table(
            file = paste0(css_dir, filename),
            sep="\t", header = TRUE)

  # Make volcano plot before filtering
  volcano(df, plot_out_dir,LM22_type,"_all_junctions")

  # Filter to most significant junction per cluster
  df_filtered <- filter_top_junction(df)
  print(dim(df))
  print(dim(df_filtered))

  write.table(df_filtered, file=paste0(nchar(filename)-4,"_top_junc_per_cluster.tsv"),sep = "\t")

  # Make volcano plot after filtering
  volcano(df_filtered, plot_out_dir,LM22_type,"_filtered_junctions")

  # Get top negative delta events
  top_sig_by_pval_negdelta <- df_filtered %>%
      dplyr::filter(p.value <= .01 ) %>%
      dplyr::filter(delta < -.2 ) %>%
      dplyr::arrange(p.value) %>%
      head(topN) %>%
      pull(event)

  # print("top_sig_by_pval_negdelta:")
  # print(top_sig_by_pval_negdelta)

  # Make plots for top negative events
  lapply(top_sig_by_pval_negdelta,  plot_event, cell_type = LM22_type,
        LM_type=meta_col, out_dir = plot_out_dir)

  # Get top positive delta events
  top_sig_by_pval_posdelta <- df_filtered %>%
    dplyr::filter(p.value <= .01) %>%
    dplyr::filter(delta > .2 ) %>%
    dplyr::arrange(p.value) %>%
    head(topN) %>%
    pull(event)


  # print("top_sig_by_pval_posdelta:")
  # print(top_sig_by_pval_posdelta)

  # Make plots for top positive events
  lapply(top_sig_by_pval_posdelta,  plot_event, cell_type = LM22_type,
          LM_type=meta_col,
          out_dir = plot_out_dir )

  return(list(
              "top_pos" = top_sig_by_pval_posdelta,
              "top_neg" = top_sig_by_pval_negdelta,
              "top_neg_and_pos" = unlist(list(top_sig_by_pval_negdelta,top_sig_by_pval_posdelta)))
              )

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
      # geom_text(stat='count', aes(label=after_stat(count)), nudge_y = -2)


  ggsave(plot = p, filename = paste0(out_dir,cell_type,"/",sig_event,".png"))

}

unpack_import_css_res <- function(ls_res){
  #' Unpack the results from the import_mesa_css() into lists
  #'
  #' @param ls_res output from
  #' @param list of top positive evnet, top negative events, top negative and positive

  ls_top_pos <- ls_top_neg <- ls_top_neg_and_pos  <- c()
  for (item in ls_res) {
       ls_top_pos <- append(ls_top_pos, item[1])
       ls_top_neg <- append(ls_top_neg, item[2])
       ls_top_neg_and_pos <- append(ls_top_neg_and_pos, item[3])
     }

  ls_top_pos<- unlist(ls_top_pos)
  ls_top_neg<- unlist(ls_top_neg)
  ls_top_neg_and_pos<- unlist(ls_top_neg_and_pos)

  return(list(ls_top_pos,ls_top_neg,ls_top_neg_and_pos))
  }


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

      print(rowname_on_off)
      print(get(rowname_on_off))
      print(paste0(rowname_on_off))

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
        paste0(opt$out_dir,"/ref_matrix/",label,"_",val,"_rowname",rowname_on_off,".pdf"))
      }
}

import_mesa_to_heatmap<- function(ls_cell_types, top_n,  label, css_dir, meta_col){
  #' Import results from MESA compare_sample_sets runs within a broader cell type
  #' (e.g. within T-cells). Find the top signficant events, make event level
  #' plots, make heatmaps of the events. This function calls import_mesa_css(),
  #' unpack_import_css_res(), and make_pheatmap()
  #'
  #' @param ls_cell_types - list of LM22 cell-types that were compared in
  #' compareWithinType.R script
  #' @param top_n - integer; how many of the top splicing events to use
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
                    topN=top_n,
                    import_mesa_css, meta_col =meta_col,
                    plot_out_dir =  paste0(opt$out_dir,"/ref_matrix/within_type/"),
                    css_dir =  css_dir)


  # Unpack top events into lists
  ls_top_events <- unpack_import_css_res(ls_res)

  # Filter metadata
  df_metadata_subset <- metadata %>%
    dplyr::filter(LM22 %in% ls_cell_types) %>%
    droplevels(.) %>%
    dplyr::arrange(Run)

  # Filter all PS
  df_all_PS <- all_PS %>%
    dplyr::select(as.vector(unlist(df_metadata_subset$Run)))

  # Make heatmap with this cell types events only within this cell types samples
  make_pheatmap(ls_top_events[[3]], paste0(label, "_diff_splicing_heatmap"),
          df_metadata_subset, df_all_PS )

  # Make heatmap with this cell types events and ALL samples
  make_pheatmap(ls_top_events[[3]], paste0(label,"_diff_splicing_heatmap_all_samples"),
          metadata, all_PS )

  return(ls_top_events)
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
# Remove rows where second column is empty (no ME Junction) and format
df_clusters_filter <- df_clusters %>% dplyr::filter(V2 != "")
df_clusters_filter$V2 <-  strsplit(as.character(df_clusters_filter$V2), ",")

# Make list of clusters by combining first col and the list in the second col
ls_clusters <- list()
for (row in 1:nrow(df_clusters_filter)) {
    event_main <- df_clusters_filter[row, "V1"]
    event_others  <- df_clusters_filter[row, "V2"]
    row_clusters <- list(unlist(append(as.vector(event_others),
                                      as.character(event_main))))
    ls_clusters <- append(ls_clusters, row_clusters )
}

# Reduce to unique clusters; remove A:B , B:A
# To DO - doesnt work ls_clusters <- unique(ls_clusters)

print(length(ls_clusters))

# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/LM22"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/LM22"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/LM6"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/LM6"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/within_type"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/within_type"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/UMAPs"))){
  dir.create(paste0(opt$out_dir,"/ref_matrix/UMAPs"),
   recursive = TRUE, showWarnings = TRUE)
}


########################################
# Import LM22 1 vs all comparisons
#######################################
#
# # Get all outputs from compare sample sets 1 vs all comparisons
# ls_lm22_css_file_names <- list.files(
#                                   paste0(opt$out_dir,
#                                   "/compare_LM22/mesa_css_outputs/"),
#                                   pattern = ".tsv")
# ls_lm22_css_file_names <- ls_lm22_css_file_names[!ls_lm22_css_file_names %in% c("heatmaps")]
#
# # Import, find signficant events and plot each one
# ls_lm22_res <- foreach(i=ls_lm22_css_file_names, .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
#     import_mesa_css(
#       filename = i,
#       topN = 10,
#       meta_col="LM22",
#       plot_out_dir = paste0(opt$out_dir,"/ref_matrix/LM22/"),
#       css_dir=paste0(opt$out_dir,"/compare_LM22/mesa_css_outputs/")
#     )
#     }
#
# # Unpack top events into lists
# ls_lm22_top_events <- unpack_import_css_res(ls_lm22_res)
# # print(ls_lm22_top_events)
#
# # Make heatmap using top events
# make_pheatmap(ls_lm22_top_events[[3]], "LM22_diff_splicing_heatmap", metadata, all_PS )

########################################
# Import LM6 1 vs all comparisons
#######################################
# Get all outputs from compare sample sets 1 vs all comparisons
ls_lm6_css_file_names <- list.files(
                        paste0(opt$out_dir,
                        "/compare_LM6/mesa_css_outputs/"),
                      pattern = ".tsv")[1]
ls_lm6_css_file_names <- ls_lm6_css_file_names[!ls_lm6_css_file_names %in% c("heatmaps")]

# Import, find signficant events and plot each one
ls_lm6_res <- foreach(i=ls_lm6_css_file_names, .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_mesa_css(
      filename = i,
      topN = 10,
      meta_col="LM22",
      plot_out_dir = paste0(opt$out_dir,"/ref_matrix/LM6/"),
      css_dir=paste0(opt$out_dir,"/compare_LM6/mesa_css_outputs/")
    )}


#
# # Unpack top events into lists
# ls_lm6_top_events <- unpack_import_css_res(ls_lm6_res)
# # print(ls_lm6_top_events)
#
# # Make heatmap using top events
# make_pheatmap(ls_lm6_top_events[[3]], "LM6_diff_splicing_heatmap", metadata, all_PS )

quit()

#######################################
# Import within cell type comparisons
######################################
T_cell_types <- list(
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory  activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta")
mon_mac_cell_types <- list(
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2")
B_cell_types <- list(
  "B cells naive",
  "B cells memory")
dendritic_cell_types <- list(
  "Dendritic cells resting",
  "Dendritic cells activated")
# mast_cell_types <- list(
#   "Mast cells resting",
#   "Mast cells activated")
# NK_cell_types <- list(
#   "NK cells resting",
#   "NK cells activated")

ls_within_cell_types <- list(
  "T_cell_types" = T_cell_types,
  "mon_mac_cell_types" = mon_mac_cell_types,
  "B_cell_types" = B_cell_types,
  "dendritic_cell_types" = dendritic_cell_types)

ls_within_res <- foreach(i=names(ls_within_cell_types),
                        .packages=c('magrittr','dplyr','ggplot2','pheatmap')
                        )%dopar%{

    import_mesa_to_heatmap(
      ls_cell_types = ls_within_cell_types[[i]],
      top_n = 10,
      label=i,
      css_dir=paste0(opt$out_dir,"/compare_within_type/mesa_css_outputs/"),
      meta_col = "LM22")
    }

names(ls_within_res) <- names(ls_within_cell_types)
print(ls_within_res[1])
print(length(ls_within_res))


##############################################
# Make reference matrices combining analysis
#############################################

# LM6 and LM22 heatmap
ls_lm6_lm22_top_events <- unique(unlist(append(
    as.list(droplevels(unname(ls_lm6_top_events[[3]]))),
    as.list(droplevels(unname(ls_lm22_top_events[[3]])))
  )))

print("ls_lm6_lm22_top_events")
# print(head(ls_lm6_lm22_top_events))
print(length(ls_lm6_lm22_top_events))

# LM6 and LM22 heatmap
make_pheatmap(ls_lm6_lm22_top_events, "LM6_and_LM22_diff_splicing_heatmap", metadata, all_PS )

# LM6 and LM22 intersection heatmap
ls_lm6_lm22_intersect <- ls_lm6_lm22_top_events <- unlist(intersect(
    as.list(droplevels(unname(ls_lm6_top_events[[3]]))),
    as.list(droplevels(unname(ls_lm22_top_events[[3]])))
  ))
print("ls_lm6_lm22_intersect")
# print(head(ls_lm6_lm22_intersect))
print(length(ls_lm6_lm22_intersect))
make_pheatmap(ls_lm6_lm22_intersect, "LM6_and_LM22_intersection_diff_splicing_heatmap", metadata, all_PS )

# Combined within cell type comparisons heatma
ls_within_cell_top_events <-  c()
for (item in ls_within_res) {
    # ls_top_pos<- append(, item[1])
    # ls_top_neg <- append(, item[2])
    ls_within_cell_top_events <- append(ls_within_cell_top_events, item[3])
   }

ls_all_within_cell_top_events <- droplevels(unname(unlist(lapply(ls_within_cell_top_events, unlist))))
print("ls_all_within_cell_top_events")
print(length(ls_all_within_cell_top_events))
make_pheatmap(ls_all_within_cell_top_events, "combined_within_cell_type_diff_splicing_heatmap", metadata, all_PS )

# LM6 and LM22 and within cell type comparisons heatmap
ls_lm6_lm22_withincelltype_top_events <- c(as.character(ls_lm6_lm22_top_events), as.character(ls_all_within_cell_top_events) )

print("ls_lm6_lm22_withincelltype_top_events")
# print(head(ls_lm6_lm22_withincelltype_top_events))
print(length(ls_lm6_lm22_withincelltype_top_events))
make_pheatmap(ls_lm6_lm22_withincelltype_top_events, "LM6_LM22_and_combined_within_cell_type_diff_splicing_heatmap", metadata, all_PS )


######################################
# UMAPs of PS using top events
######################################
# Read in merged allPS file

print(dim(all_PS))
all_PS_top_junctioins <- all_PS %>%
  dplyr::filter(rownames(all_PS) %in% ls_lm6_lm22_withincelltype_top_events)

print(dim(all_PS_top_junctioins))

# Drop genes with low variance.
allPS_var <- apply(all_PS_top_junctioins[, -1], 1, var)
print(length(allPS_var))

# For gene I used median (50% quantile) as the cutoff
# For splicing using 75% due to there being many more rows)
# trying 25
PS_param <- quantile(allPS_var, c(.25), na.rm=T)
print(PS_param)

df_allPS_filt <- all_PS_top_junctioins[allPS_var > PS_param & !is.na(allPS_var), ]
print(dim(df_allPS_filt))

# Transpose and format
df_allPS_filt_t <- as.data.frame(t(df_allPS_filt))
rownames(df_allPS_filt_t) <- colnames(df_allPS_filt)

print(dim(df_allPS_filt))
print(dim(df_allPS_filt_t))
# PCA.
all.ps.prcomp.out = as.data.frame(prcomp(na.omit(df_allPS_filt_t), center=T,  scale = T)$x)

# Making variations of UMAPs with different numbers of neighbors
lapply(c(15,20,25,30), make_umap, meta_col="data_source",
  df_PCA = all.ps.prcomp.out, out_path = "ref_matrix/UMAPs/UMAP")
lapply(c(15,20,25,30), make_umap, meta_col="LM22",
  df_PCA = all.ps.prcomp.out, out_path = "ref_matrix/UMAPs/UMAP")
lapply(c(15,20,25,30), make_umap, meta_col="sigil_general",
  df_PCA = all.ps.prcomp.out, out_path = "ref_matrix/UMAPs/UMAP")
lapply(c(15,20,25,30), make_umap, meta_col="LM6",
  df_PCA = all.ps.prcomp.out, out_path = "ref_matrix/UMAPs/UMAP")
