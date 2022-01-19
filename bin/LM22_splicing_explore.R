#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)

#############################
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

import_mesa_css <- function(filename, topN){

  # Filename to string
  LM22_type <-  substr(filename, 1, nchar(filename)-4)
  print(LM22_type
  )
  # Make output directories
  if (!dir.exists(paste0(opt$out_dir,"/explore/",LM22_type))){
    dir.create(paste0(opt$out_dir,"/explore/",LM22_type),
     recursive = TRUE, showWarnings = TRUE)
  }

  # Open mesa css
  df <- read.table(
            file = paste0(
              opt$out_dir,
              "/mesa_compare_outputs/mesa_css_outputs/",
              filename),
            sep="\t", header = TRUE)

  # Get top events by pvalue
  top_sig_by_pval <- df %>%
    dplyr::arrange(p.value) %>%
    head(2*topN)

  # Get top negative delta events
  top_sig_by_pval_negdelta <- top_sig_by_pval %>%
    dplyr::filter(delta < -.2 ) %>%
    dplyr::arrange(delta) %>%
    head(topN)%>%
    pull(event)

  # Make plots for top negative events
  lapply(top_sig_by_pval_negdelta,  plot_event, cell_type = LM22_type)

  # Get top positive delta events
  top_sig_by_pval_posdelta <- top_sig_by_pval %>%
    dplyr::filter(delta > .2 ) %>%
    dplyr::arrange(desc(delta)) %>%
    head(topN)%>%
    pull(event)

  # Make plots for top positive events
  lapply(top_sig_by_pval_posdelta,  plot_event, cell_type = LM22_type)

  return(list("top_pos" = top_sig_by_pval_posdelta,
              "top_neg" = top_sig_by_pval_negdelta,
              "top_neg_and_pos" = unlist(list(top_sig_by_pval_negdelta,top_sig_by_pval_posdelta)))
        )

}

plot_event <- function(sig_event, cell_type){

  df <- all_PS_meta %>%
    tibble::rownames_to_column(var = "event")%>%
    dplyr::filter(event %in% list(paste0(sig_event),"LM22"))%>%
    t() %>%
    as.data.frame()

  df_ <- df # copy df
  colnames(df_) <- c( "PS","LM22") #add column names from first row

  df_ <- df_[-1,] %>% # drop first row
          dplyr::filter(PS != "NaN")  #drop samples with Nan

  df_$PS <- as.numeric(df_$PS)

  p <- ggplot( df_, aes(x = LM22, y = PS, color=LM22)) +
      # geom_violin() +
      geom_jitter(position=position_jitter(0.2), alpha = 0.5) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(legend.position = "None") +
      scale_y_continuous(limits = c(0, 100))

  ggsave(plot = p, filename = paste0(
                opt$out_dir,"/explore/",cell_type,"/",sig_event,".png"))

}


# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--mesa_PS"),
    type = "character",
    default = NULL,
    help = " `mesa_allPS.tsv` input file "),

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
  dplyr::select(Run,LM22,sigil_general) %>%
  tibble::column_to_rownames("Run") %>%
  t()

all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE)
# all_PS <- tibble::rownames_to_column(all_PS, var = "event")

# add metadata to PS
all_PS_meta <- rbind(all_PS, df_sample_annotations)

# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/explore/"))){
  dir.create(paste0(opt$out_dir,"/explore/"),
   recursive = TRUE, showWarnings = TRUE)
}

# Get all outputs from compare sample sets 1 vs all comparisons
ls_css_file_names <- list.files(
                        paste0(opt$out_dir,
                        "/mesa_compare_outputs/mesa_css_outputs/"),
                      pattern = ".tsv")

ls_css_file_names <- ls_css_file_names[!ls_css_file_names %in% c("heatmaps")]
print(ls_css_file_names)

# Import ,  find signficant events , plot
ls_res <- lapply(ls_css_file_names, topN=10,  import_mesa_css)

ls_top_pos <- ls_top_neg <- ls_top_neg_and_pos  <- c()
for (item in ls_res) {
     ls_top_pos <- append(ls_top_pos, item[1])
     ls_top_neg <- append(ls_top_neg, item[2])
     ls_top_neg_and_pos <- append(ls_top_neg_and_pos, item[3])
   }

ls_top_pos<- unlist(ls_top_pos)
ls_top_neg<- unlist(ls_top_neg)
ls_top_neg_and_pos<- unlist(ls_top_neg_and_pos)


# make heatmap / reference matrix
# Filter MESA all PS file to events of interest
df_all_PS_sig_events <- all_PS %>%
  tibble::rownames_to_column('event') %>%
  dplyr::filter(event %in% ls_top_neg_and_pos) %>%
  tibble::column_to_rownames('event')

# DF to label samples(columns) with labels
df_sample_annotations <- metadata %>%
  dplyr::select(Run, LM22) %>%
  tibble::column_to_rownames("Run")

heatmap_res <- pheatmap(
  main = paste0(" "),
  df_all_PS_sig_events,
  scale = "row",
  show_rownames=F,
  show_colnames=F,
  na_col = "grey",
  annotation_col = df_sample_annotations)

save_pheatmap_pdf(
  heatmap_res,
  paste0(opt$out_dir,
        "/diff_splicing_heatmap.pdf"))
