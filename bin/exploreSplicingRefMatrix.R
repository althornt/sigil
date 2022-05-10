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
  chr = paste0("chr",strsplit(junction, ":")[[1]][1])
  range = strsplit(junction, ":")[[1]][2]
  range_start = strsplit(range, "-")[[1]][1]
  range_end = strsplit(range, "-")[[1]][2]

  # print(chr)
  # print(range_start)
  # print(range_end)

  bed_str <- paste0(chr,"\t",range_start,"\t",range_end)
  return(bed_str)
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

make_pheatmap <- function(ls_events, label, df_meta, df_PS){
  #' Make heatmap using the pheatmap package using the given events
  #' and data
  #'
  #' @param ls_events - list of events of interest to use in heatmap
  #' @param label - label to use in output file path
  #' @param df_meta - df of metadata with Run, val, and data_source, group_label, main_label columns
  #' @param df_PS - df of MESA all PS file

  # Filter MESA all PS file to events of interest
  df_all_PS_sig_events <- df_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_events) %>%
    dplyr::select(noquote(order(colnames(.)))) %>%
    tibble::column_to_rownames('event')

  # print(df_all_PS_sig_events)

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
calcGroupMed <- function(df_PS, label_type){
  str_label <- as.character(label_type)

  #transpose df, summarize to median of each group
  df_t_med <- df_all_PS_ref_events %>%
    rownames_to_column("rowname") %>%   #transpose
    gather(var, value, -rowname) %>% 
    spread(rowname, value) %>%
    filter(get(label_type) != "") %>%     #remove samples without a label grouping
    select(c(paste0(str_label),"var",df_spliceRefMatrix$event )) %>% #keep label col and junctions
    mutate_at(df_spliceRefMatrix$event, as.numeric) %>% # convert to numeric
    group_by_at(label_type) %>%
    summarise_at(vars(df_spliceRefMatrix$event), funs(median(., na.rm=TRUE))) %>%
    column_to_rownames(paste(str_label)) %>%
    as.data.frame(.) %>%
    select_if(~ !any(is.na(.))) #remove NA

  df_med <-  df_t_med  %>%
    rownames_to_column("rowname") %>% 
    gather(var, value, -rowname)  %>% 
    spread(rowname, value) %>% 
    as.data.frame(.) %>%
    column_to_rownames("var")
    
  #Remove junctions with 0 variance which break scaling
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
          paste0(opt$out_dir,"/",paste(label_type),"_med_heatmap.pdf")
          )

  return(df_med)
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
# df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)

print(head(df_all_PS))
# print(str(df_all_PS))
# quit()

metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Add metadata to splice ref matrix
df_all_PS_meta <- rbind(df_all_PS, df_sample_annotations)

df_all_PS_ref_events <- df_all_PS_meta %>%
  rownames_to_column("col") %>%
  filter((col %in% df_spliceRefMatrix$event)|(col %in% rownames(df_sample_annotations))) %>%
  column_to_rownames("col")

print(dim(df_all_PS_ref_events))

########################################################
# Make bed 
########################################################
print("Number of unique events:")
print(length(unique(df_spliceRefMatrix$event)))
ls_bed <- as.vector(sapply( unique(df_spliceRefMatrix$event), junc2bed))
writeLines(ls_bed, paste0(opt$out_dir,"/RefMatrix.bed"))

######################################################
# RefMat with medians + z-score
######################################################
# Calculate group medians and make heatmap then z-score
df_group_label_med <- calcGroupMed(df_all_PS_ref_events, "group_label")
df_group_label_med_z <- t(scale(t(df_group_label_med)))

df_main_label_med <- calcGroupMed(df_all_PS_ref_events, "main_label")
df_main_label_med_z <- t(scale(t(df_main_label_med)))

print(dim(df_group_label_med))
print(dim(df_group_label_med_z))

print(dim(df_main_label_med))
print(dim(df_main_label_med_z))

write.csv(df_group_label_med, paste0(opt$out_dir,"/group_label_med.csv") )
write.csv(df_group_label_med_z, paste0(opt$out_dir,"/group_label_med_zscore.csv") )
write.csv(df_main_label_med, paste0(opt$out_dir,"/main_label_med.csv") )
write.csv(df_main_label_med_z, paste0(opt$out_dir,"/main_label_med_zscore.csv") )

# ########################################################
# # Make subsets from gene lists
# ########################################################

# ls_cd_junctions <- df_spliceRefMatrix$event[grep("CD", df_spliceRefMatrix$overlapping)]
# print(ls_cd_junctions)

# ls_IL_junctions <- df_spliceRefMatrix$overlapping[grep("IL", df_spliceRefMatrix$overlapping)]
# print(ls_IL_junctions)

# print("List of genes:")
# print(sort(unique(df_spliceRefMatrix$overlapping)))
# df_spliceRefMatrix %>%
#   filter(overlapping == "DICER1") %>%
#   select(event)
  
# ls_immune_genes <- c("CD83","CD300A","CD44","IL32","IL7R")
# ls_immune_junctions <- df_spliceRefMatrix$event[df_spliceRefMatrix$overlapping %in% ls_immune_genes]

# print(ls_immune_junctions)

# df_LM22_med_immune_junctions <- df_LM22_med %>%
#   rownames_to_column("col") %>%
#   filter(col %in% ls_immune_junctions) %>%
#   column_to_rownames("col")

# print(df_LM22_med_immune_junctions)

# heatmap_res <- pheatmap(
#         main = paste0(" "),
#         df_LM22_med_immune_junctions,
#         scale = "row",
#         show_rownames=T,
#         show_colnames=T,
#         na_col = "grey"
#         )

# save_pheatmap_pdf(
#         heatmap_res,
#         paste0(opt$out_dir,"/CD_IL_genes_heatmap_med.pdf")
#         )

# print(df_spliceRefMatrix[df_spliceRefMatrix$overlapping %in% ls_immune_genes,])

# print(as.data.frame(df_spliceRefMatrix[df_spliceRefMatrix$event %in% ls_immune_junctions,]))





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
  inner_join(event_2_celltypes, by="event") %>%
  select(event, n)
# write.csv(df_countbygene, paste0(opt$out_dir,"/countbygene.csv"))

print(df_countbyevent)

# Plot distribution of counts per gene 
df_hist_event <- df_spliceRefMatrix %>% 
  count(event, sort = TRUE) 
ggplot(df_hist_event, aes(x=n)) + 
  geom_histogram() + 
  scale_x_continuous(breaks = round(seq(min(df_hist_event$n), max(df_hist_event$n), by = 1),1)) +
  ggtitle("Counts per event in the splicing reference matrix") 
ggsave( paste0(opt$out_dir,"/hist_count_per_event.png"))





