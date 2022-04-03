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
library(UpSetR)

cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)


########################
# Functions
#########################
format_merge <- function(df_gene,  df_splice){
  # Making long dfs
  df_gene_z_long <- df_gene %>%
      tidyr::gather(., key="cell_ID", value = gene_z, -c(gene, cell_types)) %>%
      as.data.frame()

  df_splice_z_long <- df_splice %>%
      tidyr::gather(., key="cell_ID", value = splice_z, -c(gene, event, cell_types)) %>%
      as.data.frame()

  # Merging the long dfs 
  df_merged_z_long <-  df_splice_z_long %>%
          inner_join(df_gene_z_long, by = c("gene","cell_ID"), suffix=c("_splice","_gene")) %>% 
          tidyr::drop_na()

  df_merged_z_long$splice_z = as.numeric(as.character(df_merged_z_long$splice_z))
  df_merged_z_long$gene_z = as.numeric(as.character(df_merged_z_long$gene_z))

  return(df_merged_z_long)

}
scatter_per_gene <- function(df_long, ouput_prefix){

  # Make scatter plot for each matching gene
  ls_genes <- unique(df_long$gene)
  for (g in ls_genes){
    df_g <- df_long %>%
      filter(gene == g)
    p <- ggplot(aes(x=gene_z, y=splice_z, color = cell_ID, shape = event), data=df_g)+ 
      geom_point() +
      labs(title= g) + geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
      geom_text(
                label= df_g$cell_types_splice,
                nudge_x = 0.05, nudge_y = 0.05,
                check_overlap =F, col = "black", size = 1
              ) +
      theme_minimal()

    ggsave(plot = p, filename = paste0(ouput_prefix, g, ".png"))
  }

}
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#   #' Function to save pheatmaps to a pdf file
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }

scatter_all <- function(df_long, output_path){
  # Scatter plot with all matching events/genes
  p <- ggplot(aes(x=gene_z, y=splice_z, color = gene), data=df_long)+ 
    geom_point() +
    geom_text(
              label= df_long$gene,
              nudge_x = 0.05, nudge_y = 0.05,
              check_overlap =F, col = "black", size = 1.5
            )
  ggsave(plot = p, filename = output_path)
}

calcGroupMed <- function(df_exp, ls_gene, LM_type){
  str_LM <- as.character(LM_type)
  meta_row <- df_sample_annotations[c(paste0(LM_type)),]

  # Add LM type metadata row to expression
  df_exp_meta <- rbind(df_exp, meta_row )
  rownames(df_exp_meta)[length(rownames(df_exp_meta))] <-paste0(LM_type) # name new row
  print(ls_gene)

  ls_gene <- intersect(ls_gene, rownames(df_exp))

  print(ls_gene)

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
    
  # #Remove genes with 0 variance which break scaling
  # df_med_var<- df_med[apply(df_med, 1, var) != 0, ]

  # heatmap_res <- pheatmap(
  #         main = paste0(" "),
  #         df_med_var,
  #         scale = "row",
  #         show_rownames=F,
  #         show_colnames=T,
  #         na_col = "grey"
  #         )

  # save_pheatmap_pdf(
  #         heatmap_res,
  #         paste0(opt$out_dir,"/",paste(LM_type),"_med_heatmap.pdf")
  #         )

  return(df_med)
}

# Arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--spliceDir"),
    type = "character",
    default = NULL,
    help = " "),

  optparse::make_option(
    c("-c", "--geneDir"),
    type = "character",
    default = NULL,
    help = " "),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"), 
  
  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to put outputs")
    )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

print(opt$spliceDir)
print(opt$geneDir)

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/LM6_scatterplots/"))){
  dir.create(file.path(opt$out_dir,"/LM6_scatterplots/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/LM22_scatterplots/"))){
  dir.create(file.path(opt$out_dir,"/LM22_scatterplots/"),
              recursive = TRUE, showWarnings = TRUE)}

# Read in files 

# Metadata
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Events in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 
print("splice ref mat")
print(head(df_splice_ref))
print(dim(df_splice_ref))

print(unique(df_splice_ref$group))

# Map events to their signifcant cell type
df_event2cell <- df_splice_ref %>%
    group_by(event) %>%
    summarize(context = list(cell_type)) %>%
    mutate(cell_types = map_chr(context, toString)) %>%
    select(event, cell_types) %>%
    as.data.frame()
print(head(df_event2cell))
print(dim(df_event2cell))

df_event2cell <- merge(x = df_event2cell, y = df_splice_ref, by="event") %>%
  distinct(event, .keep_all = TRUE)

print(head(df_event2cell))
print(dim(df_event2cell))

print("=====================================================================")

# Genes in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(
                              opt$geneDir,
                              "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) %>%
              rename(gene = X)
print("gene ref mat")
print(head(df_gene_ref))

print(unique(df_gene_ref$group))
print(dim(df_gene_ref))


# Map genes to to their signifcant cell type
df_DEG2cell <- df_gene_ref %>%
    group_by(gene) %>%
    summarize(context = list(cell_type)) %>%
    mutate(cell_types = map_chr(context, toString)) %>%
    select(gene, cell_types) %>%
    as.data.frame()

print(head(df_DEG2cell))
print(dim(df_DEG2cell))

#######################
# LM6 group z_scores
########################

# Open LM6 splice z scores 
df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix/LM6_med_zscore.csv"),
                            header=TRUE) %>%
                    rename(event = X)

# Add sig cell type and gene name
df_lm6_splice_z <-  merge(x=df_lm6_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
    select( -column_label2,-group, -column_label, -cell_type) %>%
    rename(gene = overlapping)

print("----------------------------------------")
print(head(df_lm6_splice_z))
print(dim(df_lm6_splice_z))

# Open LM6 gene z scores 
df_lm6_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE) %>%
    rename(gene = X)
# Add  cell_type to z scores
df_lm6_gene_z <- merge(x=df_DEG2cell,y=df_lm6_gene_z,by="gene",all=TRUE) 
print("df_lm6_gene_z")
print(head(df_lm6_gene_z))
print(dim(df_lm6_gene_z))


# # Scatterplots 
# df_lm6_merged_z_long <- format_merge(df_lm6_gene_z,df_lm6_splice_z )
# print(head(df_lm6_merged_z_long))
# scatter_all(df_lm6_merged_z_long, paste0(opt$out_dir, "/lm6_scatter_all.png"))
# scatter_per_gene(df_lm6_merged_z_long, paste0(opt$out_dir, "/LM6_scatterplots/"))

# #######################
# # LM22
# ########################
# Open LM22 splice z scores 
df_lm22_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix/LM22_med_zscore.csv"),
                            header=TRUE) %>%
                    rename(event = X)

# Add sig cell type and gene name
df_lm22_splice_z <-  merge(x=df_lm22_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
    select( -column_label2,-group, -column_label, -cell_type) %>%
    rename(gene = overlapping)

# print("----------------------------------------")
# print(head(df_lm22_splice_z))
# print(dim(df_lm22_splice_z))

# # Open LM22 gene z scores 
# df_lm22_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM22_med_zscore.csv"), header=TRUE) %>%
#     rename(gene = X)
# # Add  cell_type to z scores
# df_lm22_gene_z <- merge(x=df_DEG2cell,y=df_lm22_gene_z,by="gene",all=TRUE) 
# print("df_lm22_gene_z")
# print(head(df_lm22_gene_z))
# print(dim(df_lm22_gene_z))

# # # Scatterplots 
# df_lm22_merged_z_long <- format_merge(df_lm22_gene_z,df_lm22_splice_z )
# scatter_all(df_lm22_merged_z_long, paste0(opt$out_dir, "/lm22_scatter_all.png"))
# scatter_per_gene(df_lm22_merged_z_long, paste0(opt$out_dir, "/LM22_scatterplots/"))

##############################
# Count splicing types 
##############################
	# - event - NA - no matched gene
	# - event - matched gene not in gene ref matrix
	# - event - matched gene in ref matrix same cell 
	# - event - matched gene in ref matrix different cell 
# print("Total number of events:")
# print(nrow(df_splice_ref))

# print("Total number of unique events:")
# print(length(unique(df_splice_ref$event)))

# print("Number of events with no gene:")
# print( nrow(df_splice_ref %>% 
#   filter(overlapping == "")))

# df_merged_z <- df_lm22_splice_z  %>%
#         inner_join(df_lm22_gene_z, by = c("gene"), suffix=c("_splice","_gene")) %>%
#         select(event, gene, cell_types_splice, gene,cell_types_gene ) %>%
#         arrange(gene)
        
# print("Number of events with matched gene in gene reference matrix:")
# print(nrow(df_merged_z))

# print("Number of spliced genes with matched gene in gene reference matrix:")
# print(length(unique(df_merged_z$gene)))

##############################
# LM22 splice ref vs all matched genes
##############################
# Read in gene exp 
df_exp <- read.table(file= paste0(
                              opt$geneDir,
                              "combined_kallisto_log2tpm_batch_corrected.csv"),
                    sep=",", header = TRUE, row.names=1)
df_exp <- df_exp %>% mutate_if(is.character,as.numeric)
print(head(df_exp))
print(dim(df_exp))

# get medians and z-score across group
df_exp_LM22_med <- calcGroupMed(df_exp, unlist(unique(df_splice_ref$overlapping)), "LM22")
df_exp_LM22_med_z <- t(scale(t(df_exp_LM22_med)))
df_exp_LM6_med <- calcGroupMed(df_exp, unlist(unique(df_splice_ref$overlapping)), "LM6")
df_exp_LM6_med_z <- t(scale(t(df_exp_LM6_med)))

print(head(df_exp_LM22_med_z))
print(dim(df_exp_LM22_med_z))
print("------")
print(head(df_exp_LM6_med_z))
print(dim(df_exp_LM6_med_z))

# Combine LM6 and LM22
df_exp_LM6_LM22_med_z <- cbind(df_exp_LM6_med_z, df_exp_LM22_med_z)
print("Combined......")
print(head(df_exp_LM6_LM22_med_z))
print(dim(df_exp_LM6_LM22_med_z))

# Replace spaces with _ to match other df
str_conv <- function(in_str){
  paste(unlist(strsplit(as.character(in_str), split=" ")),
                                    collapse="_")
}

# Replace spaces with _ to match other df
new_names <- unlist(lapply(colnames(df_exp_LM6_LM22_med_z),FUN= str_conv))
colnames(df_exp_LM6_LM22_med_z) <- new_names
print(head(df_exp_LM6_LM22_med_z))

# Check new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type %in% colnames(df_exp_LM6_LM22_med_z))))

# Go through each splice event and get row for that gene 
print(head(df_splice_ref))

# # All genes that were in splice ref matrix
df_splice_ref_genes <- unique(df_splice_ref$overlapping)
print(length(df_splice_ref_genes))
print("=======================================================================")

df_splice_ref['gene_z'] <- NA
df_splice_ref['splice_z'] <- NA
print(head(df_splice_ref))
# quit()
# Adding gene z_scores to splice ref df
# for (gene in head(df_splice_ref_genes)){

for (gene in df_splice_ref_genes){
  print("start.............")
  print(gene)
  #  For each gene get the cell_types the splice event was signifcant in
  df_  <- df_splice_ref %>% filter(overlapping == gene)
  # print(as.vector(unlist(unique(df_$cell_type))))

  # Check if gene in gene exp data (splicing can have lists for overlaps "AC129492.1,AC129492.4")
  if (gene %in% rownames(df_exp_LM6_LM22_med_z)){
    
    # Iterate over cell types
    for (gene_cell in as.vector(unlist(unique(df_$cell_type)))){
      # print(gene_cell)

      # Look up the gene exp val for gene and cell type
      z <- df_exp_LM6_LM22_med_z[gene, gene_cell] 
      # print(z)

      # Add gene z to each corresponding splice ref rows
      df_splice_ref <- df_splice_ref %>%
        mutate(gene_z = ifelse(((overlapping==gene) & (cell_type==gene_cell)), z, gene_z))
    }
  } else
  {print("no gene match")}
}


print(head(df_splice_ref))
# print(df_splice_ref)

# for (event in df_splice_ref_genes$event){




# }


# Adding splice z_scores to splice ref df


# event   cell type gene_zscore splice_zscore 



	# - for each reference matrix splice junction
	# 	- for that cell type pull the gene expression from the big table
	# 	- z-score all
	# 	- scatter plot one labeling gene
	# 	- another labeling event
	# 	- another labeling cell-type sig
	# 	- only plot the z-score for the cell_type specific one 

