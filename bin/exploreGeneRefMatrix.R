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
library(UpSetR)
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  #' Function to save pheatmaps to a pdf file
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
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
        paste0(opt$out_dir,"/",label,"_",val,"_rowname",rowname_on_off,".pdf"))
      }
}

calcGroupMed <- function(df_exp, ls_gene, LM_type){
  str_LM <- as.character(LM_type)
  meta_row <- df_sample_annotations[c(paste0(LM_type)),]

  # Add LM type metadata row to expression
  df_exp_meta <- rbind(df_exp, meta_row )
  rownames(df_exp_meta)[length(rownames(df_exp_meta))] <-paste0(LM_type) # name new row

  ls_gene <- intersect(ls_gene, rownames(df_exp))

  # transpose df, summarize to median of each group
  df_t_med <- df_exp_meta %>%
    rownames_to_column("row_name") %>%  
    filter((row_name %in% ls_gene) | (row_name %in% c(paste0(LM_type)))) %>%  
    gather(var, value, -row_name) %>% 
    spread(row_name, value) %>%
    filter(get(LM_type) != "") %>%     #remove samples without a LM grouping
    select(c(paste0(str_LM),"var", ls_gene)) %>% #keep LM col and gene subset
    mutate_at(vars(as.vector(ls_gene)), as.numeric) %>% # convert to numeric
    group_by_at(LM_type) %>%
    summarise_at(vars(ls_gene), funs(median(., na.rm=TRUE))) %>%
    column_to_rownames(paste(str_LM)) %>%
    as.data.frame(.) #%>%
  #   # select_if(~ !any(is.na(.))) #remove NA

  df_med <-  df_t_med  %>%
    rownames_to_column("rowname") %>% 
    gather(var, value, -rowname)  %>% 
    spread(rowname, value) %>% 
    as.data.frame(.) %>%
    column_to_rownames("var")
    
  #Remove genes with 0 variance which break scaling
  df_med_var<- df_med[apply(df_med, 1, var) != 0, ]

  heatmap_res <- pheatmap(
          main = paste0(" "),
          df_med_var,
          scale = "row",
          show_rownames=F,
          show_colnames=T,
          na_col = "grey"
          )

  save_pheatmap_pdf(
          heatmap_res,
          paste0(opt$out_dir,"/",paste(LM_type),"_med_heatmap.pdf")
          )

  return(df_med)
}

# Arguments
option_list <- list(
  optparse::make_option(
    c("-s", "--geneRefMatrix"),
    type = "character",
    default = NULL,
    help = "gene reference matrix"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),
      
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = " `combined and batch corrected input file ")
    )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
if (!dir.exists(paste0(opt$out_dir))){
  dir.create(paste0(opt$out_dir),
   recursive = TRUE, showWarnings = TRUE)
}

# Open files
# File where each row is top N result from a comparison 
df_geneRefMatrix <- read.table(file = opt$geneRefMatrix, header = TRUE)
print(head(df_geneRefMatrix))
print(dim(df_geneRefMatrix))

# Batch corrected gene expression 
df_exp <- read.table(file = opt$i, sep=",", header = TRUE, row.names=1) 
df_exp <- df_exp %>% mutate_if(is.character,as.numeric)
print(dim(df_exp))

# Metadata
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Add metadata to gene expression
# df_exp_meta <- rbind(df_exp, df_sample_annotations)
# print(dim(df_exp_meta))

# # Gene exp for the ref events 
# df_exp_ref_events <- df_exp_meta %>%
#   rownames_to_column("col") %>%
#   filter((col %in% df_geneRefMatrix$X)|(col %in% rownames(df_sample_annotations))) %>%
#   column_to_rownames("col")

# print(tail(df_exp_ref_events))
# print(dim(df_exp_ref_events))
######################################################
# RefMat with medians + z-score
######################################################

# Calculate group medians and make heatmap then z-score
df_LM6_med <- calcGroupMed(df_exp, unlist(df_geneRefMatrix$X), "LM6")
df_LM6_med_z <- t(scale(t(df_LM6_med)))

df_LM22_med <- calcGroupMed(df_exp, unlist(df_geneRefMatrix$X), "LM22")
df_LM22_med_z <- t(scale(t(df_LM22_med)))

print(dim(df_LM6_med))
print(dim(df_LM6_med_z))

print(dim(df_LM22_med))
print(dim(df_LM22_med_z))

write.csv(df_LM6_med, paste0(opt$out_dir,"/LM6_med.csv"))
write.csv(df_LM6_med_z, paste0(opt$out_dir,"/LM6_med_zscore.csv"))
write.csv(df_LM22_med, paste0(opt$out_dir,"/LM22_med.csv") )
write.csv(df_LM22_med_z, paste0(opt$out_dir,"/LM22_med_zscore.csv"))
######################################################
# Count gene occurances in matrix
######################################################

# print(head(df_geneRefMatrix))

# df_geneRefMatrix_count <- df_geneRefMatrix %>%
#     group_by(X) %>%
#     summarize(context = list(cell_type)) %>%
#   mutate(cell_types = map_chr(context, toString)) %>%
#   select(X,cell_types) 

# write.csv(df_geneRefMatrix_count, 
#             paste0(opt$out_dir,"/gene_occurances.csv"))

# # Plot distribution of counts per gene 
# df_hist_gene <- df_geneRefMatrix %>% 
#   count(X, sort = TRUE) 
# ggplot(df_hist_gene, aes(x=n)) + 
#   geom_histogram() + 
#   scale_x_continuous(breaks = round(seq(min(df_hist_gene$n), max(df_hist_gene$n), by = 1),1)) +
#   ggtitle("Counts per gene in the gene reference matrix") 
# ggsave( paste0(opt$out_dir,"/hist_count_per_gene.png"))

# ######################################################
# # Z-score all genes
# ######################################################

# this takes way too llong 
# Calculate group medians and make heatmap then z-score
# df_LM6_med_all_genes <- calcGroupMed(df_exp, unlist(rownames(df_exp)), "LM6")
# df_LM6_med_z_all_genes <- t(scale(t(df_LM6_med__all_genes)))

# df_LM22_med_all_genes <- calcGroupMed(df_exp, unlist(rownames(df_exp)), "LM22")
# df_LM22_med_z__all_genes <- t(scale(t(df_LM22_med__all_genes)))

# print(dim(df_LM6_med_all_genes))
# print(dim(df_LM6_med_z_all_genes))

# print(dim(df_LM22_med_all_genes))
# print(dim(df_LM22_med_z_all_genes))

# write.csv(df_LM6_med_all_genes, paste0(opt$out_dir,"/LM6_med_all_genes.csv"))
# write.csv(df_LM6_med_z_all_genes, paste0(opt$out_dir,"/LM6_med_zscore_all_genes.csv"))
# write.csv(df_LM22_med_all_genes, paste0(opt$out_dir,"/LM22_med_all_genes.csv") )
# write.csv(df_LM22_med_z_all_genes, paste0(opt$out_dir,"/LM22_med_zscore_all_genes.csv"))