#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)

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

import_mesa_css <- function(filename, topN, plot_out_dir){

  # Filename to string
  LM22_type <-  substr(filename, 1, nchar(filename)-4)
  # print(LM22_type)
  # Make output directories
  if (!dir.exists(paste0(plot_out_dir,LM22_type))){
    dir.create(paste0(plot_out_dir,LM22_type),
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
    dplyr::filter(delta < -.1 ) %>%
    dplyr::arrange(delta) %>%
    head(topN)%>%
    pull(event)

  # Make plots for top negative events
  lapply(top_sig_by_pval_negdelta,  plot_event, cell_type = LM22_type,
    out_dir = plot_out_dir)

  # Get top positive delta events
  top_sig_by_pval_posdelta <- top_sig_by_pval %>%
    dplyr::filter(delta > .1 ) %>%
    dplyr::arrange(desc(delta)) %>%
    head(topN)%>%
    pull(event)

  # Make plots for top positive events
  lapply(top_sig_by_pval_posdelta,  plot_event, cell_type = LM22_type,
          out_dir = plot_out_dir )

  return(list("top_pos" = top_sig_by_pval_posdelta,
              "top_neg" = top_sig_by_pval_negdelta,
              "top_neg_and_pos" = unlist(list(top_sig_by_pval_negdelta,top_sig_by_pval_posdelta)))
        )
}

plot_event <- function(sig_event, cell_type, out_dir){

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

  ggsave(plot = p, filename = paste0(out_dir,cell_type,"/",sig_event,".png"))

}

unpack_import_css_res <- function(ls_res){
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
all_PS_meta <- rbind(all_PS, df_sample_annotations)

# Make output directories
if (!dir.exists(paste0(opt$out_dir,"/explore/LM22"))){
  dir.create(paste0(opt$out_dir,"/explore/LM22"),
   recursive = TRUE, showWarnings = TRUE)
}

if (!dir.exists(paste0(opt$out_dir,"/explore/within_type"))){
  dir.create(paste0(opt$out_dir,"/explore/within_type"),
   recursive = TRUE, showWarnings = TRUE)
}

#########################################
# Import LM22 1 vs all comparisons
########################################
# Get all outputs from compare sample sets 1 vs all comparisons
ls_css_file_names <- list.files(
                        paste0(opt$out_dir,
                        "/mesa_compare_outputs/mesa_css_outputs/"),
                      pattern = ".tsv")
ls_css_file_names <- ls_css_file_names[!ls_css_file_names %in% c("heatmaps")]

# Import, find signficant events and plot each one
ls_res <- lapply(ls_css_file_names, topN=10,  import_mesa_css, plot_out_dir = paste0(opt$out_dir,"/explore/LM22/"))
# Unpack top events into lists
ls_top_events <- unpack_import_css_res(ls_res)

# ls_top_pos <- ls_top_neg <- ls_top_neg_and_pos  <- c()
# for (item in ls_res) {
#      ls_top_pos <- append(ls_top_pos, item[1])
#      ls_top_neg <- append(ls_top_neg, item[2])
#      ls_top_neg_and_pos <- append(ls_top_neg_and_pos, item[3])
#    }
#
# ls_top_pos<- unlist(ls_top_pos)
# ls_top_neg<- unlist(ls_top_neg)
# ls_top_neg_and_pos<- unlist(ls_top_neg_and_pos)

#############################
# heatmap all diff splicng
#############################

# make heatmap / reference matrix
# Filter MESA all PS file to events of interest
df_all_PS_sig_events <- all_PS %>%
  tibble::rownames_to_column('event') %>%
  dplyr::filter(event %in% ls_top_events[[3]]) %>%
  tibble::column_to_rownames('event')

for (val in list("LM22", "LM6"))
{
# DF to label samples(columns) with labels
df_sample_annotations <- metadata %>%
  dplyr::select(Run, val) %>%
  tibble::column_to_rownames("Run")

heatmap_res <- pheatmap(
  main = paste0(" "),
  df_all_PS_sig_events,
  # scale = "row",
  show_rownames=F,
  show_colnames=F,
  na_col = "grey",
  annotation_col = df_sample_annotations)

save_pheatmap_pdf(
  heatmap_res,
  paste0(opt$out_dir,
        "/explore/diff_splicing_heatmap_",val,".pdf"))
}

#########################################
# Import T-cell 1 vs all comparisons
########################################
# Get output files from compareWithinType script
ls_css_file_names_tcells <- list.files(
                        paste0(opt$out_dir,
                        "/compareWithinType/mesa_compare_outputs/mesa_css_outputs/"),
                        pattern = ".tsv")

# Import files, find top signficant events, plot each event
ls_res_tcells <- lapply(
                        ls_css_file_names_tcells,
                        topN=10,
                        import_mesa_css,
                        plot_out_dir =  paste0(opt$out_dir,"/explore/within_type/"))

# Unpack top events into lists
tcell_top_events <- unpack_import_css_res(ls_res_tcells)

print(tcell_top_events[[3]])

#############################
# heatmap T-cell
#############################

T_cell_types <- list(
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory  activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta")

# Get samples with this cell type
metadata_T_cell_types <- metadata %>%
  dplyr::filter(LM22 %in% T_cell_types)

# Filter MESA all PS file to events from T-cell types
df_all_PS_sig_events_tcell <- all_PS[metadata_T_cell_types$Run]
df_all_PS_sig_events_tcell <- df_all_PS_sig_events_tcell %>%
  tibble::rownames_to_column('event') %>%
  dplyr::filter(event %in% tcell_top_events[[3]]) %>%
  tibble::column_to_rownames('event')


print(colnames(df_all_PS_sig_events_tcell))

for (val in list("LM22", "LM6"))
{
  print(val)
  # DF to label samples(columns) with labels
  df_sample_annotations <- metadata_T_cell_types %>%
    dplyr::select(Run, val) %>%
    tibble::column_to_rownames("Run")

  print(df_sample_annotations)

  heatmap_res <- pheatmap(
    main = paste0(" "),
    df_all_PS_sig_events_tcell,
    # scale = "row",
    show_rownames=T,
    show_colnames=F,
    na_col = "grey",
    annotation_col = df_sample_annotations)

  save_pheatmap_pdf(
    heatmap_res,
    paste0(opt$out_dir,
          "/explore/diff_splicing_heatmap_",val,"_tcell.pdf"))
}



############################
# monocytes and macrophages
##################################
mon_mac_cell_types <- list(
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2")
