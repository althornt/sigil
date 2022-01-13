#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(dplyr)
library(pheatmap)
library(ggplot2)
# library(tidyr)

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

ls_css_outputs <- list.files(paste0(opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/"))
print(ls_css_outputs)

df_dend_rest <-  read.table(file = paste0(
                  opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/Dendritic_cells_resting.tsv"),
                  sep="\t", header = TRUE)
df_dend_act <-  read.table(file = paste0(
                  opt$out_dir, "/mesa_compare_outputs/mesa_css_outputs/Dendritic_cells_activated.tsv"),
                  sep="\t",  header = TRUE)

df_dend_rest %>%
  dplyr::arrange(corrected) %>%
  head()

str_event_example <- as.character(df_dend_act %>%
  dplyr::arrange(corrected) %>%
  head(1) %>%
  pull(event))


# make df for the 1 example event
event_ex <- all_PS_meta %>%
  tibble::rownames_to_column(var = "event") %>%
  dplyr::filter(event %in% list(str_event_example,"LM22")) %>%
  t() %>%
  as.data.frame()

new_event_ex <- event_ex
colnames(new_event_ex) <- c( "PS","LM22") #add column names from first row
new_event_ex <- new_event_ex[-1,] # drop first row

print(new_event_ex)

new_event_ex %>%
  dplyr::group_by(LM22)


# p <-    ggplot( new_event_ex, aes(PS)) +
#     geom_freqpoly(  binwidth = 5, stat = "count") +
#     labs(fill="")
#
# ggsave(plot = p, filename = paste0(opt$out_dir,"/explore/test.png"))

# plot
# r5roup?
