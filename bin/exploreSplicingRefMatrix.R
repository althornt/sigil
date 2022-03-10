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
junc2bed <- function(junction){
  chr = strsplit(junction, ":")[[1]][1]
  range = strsplit(junction, ":")[[1]][2]
  range_start = strsplit(range, "-")[[1]][1]
  range_end = strsplit(range, "-")[[1]][2]

  # print(chr)
  # print(range_start)
  # print(range_end)

  bed_str <- paste0(chr,"\t",range_start,"\t",range_end)
  return(bed_str)
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
    help = "full path to put outputs"),
      
  optparse::make_option(
    c("-i", "--mesa_PS"),
    type = "character",
    default = NULL,
    help = " `mesa_allPS.tsv` input file ")
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
df_spliceRefMatrix <- readr::read_tsv(file = opt$spliceRefMatrix)
print(head(df_spliceRefMatrix))
print(dim(df_spliceRefMatrix))

df_all_PS <- read.table(file = opt$mesa_PS, sep="\t", header = TRUE, row.names=1)
print(head(df_all_PS))


metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# # Add metadata to splice ref matrix
# df_all_PS_meta <- rbind(df_all_PS, df_sample_annotations)

# print(head(df_all_PS_meta))
# # print(rownames(df_all_PS_meta))

# df_all_PS_ref_events <- df_all_PS_meta %>%
#   rownames_to_column("col") %>%
#   filter((col %in% df_spliceRefMatrix$event)|(col %in% rownames(df_sample_annotations))) %>%
#   column_to_rownames("col")

# print(dim(df_all_PS_ref_events))
# # print(tail(df_all_PS_ref_events))


# ########################################################
# # RefMat with medians
# ####################################################
# df_med <- df_all_PS_ref_events %>%
#   rownames_to_column("rowname") %>% 
#   gather(var, value, -rowname) %>% 
#   spread(rowname, value) %>%
#   group_by(LM22)  
  
#   # %>%
#   # summarize(med_PS = median(LM22, na.rm = TRUE))
#   # summarise(A = median(A[!is.na(a)]), B = median(B[!is.na(b)]))

# head(df_med)


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
print(df_countbygene)

# Plot distribution of counts per gene 
df_hist <- df_spliceRefMatrix %>% 
  count(overlapping, sort = TRUE) 
ggplot(df_hist, aes(x=n)) + 
  geom_histogram() + 
  scale_x_continuous(breaks = round(seq(min(df_hist$n), max(df_hist$n), by = 1),1)) +
  ggtitle("Counts per gene in the splicing reference matrix") 
ggsave( paste0(opt$out_dir,"/hist_count_per_gene.png"))

# print("---------------------------------------------------------------------------------")


# Map each gene to cell types the event is present in 
event_2_celltypes <- df_spliceRefMatrix %>% 
  group_by(event) %>% 
  summarize(context = list(cell_type)) %>%
  mutate(cell_types = map_chr(context, toString)) %>%
  select(event,cell_types) 

print(event_2_celltypes)

# count genes in matrix
df_countbyevent <- df_spliceRefMatrix %>% 
  count(event, sort = TRUE) %>%
  inner_join(event_2_celltypes, by="event")
# write.csv(df_countbygene, paste0(opt$out_dir,"/countbygene.csv"))

print(df_countbyevent)

# Plot distribution of counts per gene 
df_hist_event <- df_spliceRefMatrix %>% 
  count(event, sort = TRUE) 
ggplot(df_hist_event, aes(x=n)) + 
  geom_histogram() + 
  scale_x_continuous(breaks = round(seq(min(df_hist$n), max(df_hist$n), by = 1),1)) +
  ggtitle("Counts per event in the splicing reference matrix") 
ggsave( paste0(opt$out_dir,"/hist_count_per_event.png"))

print("---------------------------------------------------------------------------------")




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



match('CD',df_spliceRefMatrix$overlapping)
# any(df_spliceRefMatrix$overlapping %like% 'CD')



ls_cd_junctions <- df_spliceRefMatrix$overlapping[grep("CD", df_spliceRefMatrix$overlapping)]
print(ls_cd_junctions)



ls_IL_junctions <- df_spliceRefMatrix$overlapping[grep("IL", df_spliceRefMatrix$overlapping)]
print(ls_IL_junctions)


print(sort(unique(df_spliceRefMatrix$overlapping)))

df_spliceRefMatrix %>%
  filter(overlapping == "DICER1") %>%
  select(event)
  
########################################################
# Make bed 
########################################################

print(unique(df_spliceRefMatrix$event))

# for (i in unique(df_spliceRefMatrix$event)){  

#   sink(paste0(opt$out_dir,"/spliceRefMatrix.bed"))                                   # Apply sink & cat functions

#   chr = strsplit(i, ":")[[1]][1]
#   range = strsplit(i, ":")[[1]][2]
#   range_start = strsplit(range, "-")[[1]][1]
#   range_end = strsplit(range, "-")[[1]][2]
  
#   cat(paste0(chr,"\t",range_start,"\t",range_end))
#   cat("\n")
#   close(paste0(opt$out_dir,"/spliceRefMatrix.bed"))                                   # Apply sink & cat functions


# }




ls_bed <- as.vector(sapply( unique(df_spliceRefMatrix$event), junc2bed))
writeLines(ls_bed, paste0(opt$out_dir,"/spliceRefMatrix.bed"))