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

################
# Functions
################
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

  # print(df_all_PS_sig_events)

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
        paste0(opt$out_dir,"/ref_matrix/",label,"_",val,"_rowname",rowname_on_off,".pdf"))
      }
}

# Arguments
option_list <- list(
  optparse::make_option(
    c("-s", "--spliceRefMatrix"),
    type = "character",
    default = NULL,
    help = "splicing reference matrix"),

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
df_spliceRefMatrix <- readr::read_tsv(file = opt$spliceRefMatrix)
print(head(df_spliceRefMatrix))
print(dim(df_spliceRefMatrix))

# metadata <- read.csv(file = opt$metadata)
# df_sample_annotations <- metadata %>%
#   dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
#   tibble::column_to_rownames("Run") %>%
#   t()

# # Add metadata to splice ref matrix
# df_spliceRefMatrix_meta <- rbind(df_spliceRefMatrix, df_sample_annotations)

# df_junc2gene <- read.table(file = opt$junc2gene, sep="\t",
#                             header = TRUE, row.names=1)
# head(df_junc2gene)
########################################################
# Count distribution of junctions for each gene
########################################################

print("Most common genes in the reference matrix:")

# Map each gene to cell types the event is present in 
gene_2_celltypes <- df_spliceRefMatrix %>% 
  group_by(overlapping) %>% 
  summarize(context = list(cell_type)) %>%
  mutate(cell_types = map_chr(context, toString)) %>%
  select(overlapping,cell_types) 

print(gene_2_celltypes)

# count genes in matrix
df_countbygene <- df_spliceRefMatrix %>% 
  count(overlapping, sort = TRUE) %>%
  inner_join(gene_2_celltypes, by="overlapping")
write.csv(df_countbygene, paste0(opt$out_dir,"/countbygene.csv"))

# Plot distribution of counts per gene 
df_hist <- df_spliceRefMatrix %>% 
  count(overlapping, sort = TRUE) 
ggplot(df_hist, aes(x=n)) + 
  geom_histogram() + 
  ggtitle("Number of events per gene in the splicing reference matrix") 
ggsave( paste0(opt$out_dir,"/hist_count_per_gene.png"))

print("---------------------------------------------------------------------------------")



print("Most common events in the reference matrix:")
df_spliceRefMatrix %>% 
  count(event,overlapping, sort = TRUE) 

# df_spliceRefMatrix %>% 
#   count(event,cell_type, sort = TRUE) 

# df_spliceRefMatrix %>% 
#   group_by(overlapping) %>% 
#   summarize(context = list(event))

# print("---------------------------------------------------------------------------------")



# df_spliceRefMatrix %>% 
#   group_by(overlapping) %>% 
#   summarize(context = list(cell_type))

# print(typeof(df_spliceRefMatrix$cell_type))

########################################################
# Make subsets from gene lists
########################################################

# Look at all "CD" genes
# # Find junctions from genes that contain CD
# ls_cd_junctions <- ls_lm22_top_events[grep("CD", names(ls_lm22_top_events))]
# print(ls_cd_junctions)
#
# if (!dir.exists(paste0(opt$out_dir,"/ref_matrix/LM22/CD_genes"))){
#   dir.create(paste0(opt$out_dir,"/ref_matrix/LM22/CD_genes"),
#    recursive = TRUE, showWarnings = TRUE)
# }
# lapply(ls_cd_junctions,  plot_event, cell_type = "",
#       LM_type="LM22", out_dir = paste0(opt$out_dir,"/ref_matrix/LM22/CD_genes/"))


# Look at "LINC"


# Look at ""


########################################################
# Make RefMatrix with median of each  LM22 type
########################################################


# Get median of each LM22 type
# Z-score of medians
# Heatmap of medians
# make_pheatmap(ls_lm22_top_events, "LM22_diff_splicing_heatmap", metadata, all_PS )


# Get median of each LM6 type
# Z-score of medians
# Heatmap of medians
