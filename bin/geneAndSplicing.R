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

make_umap <- function(num_neighbor,meta_col,df_PCA,out_path, plot_name) {

  set.seed(123)

  # Make color palette
  n <- length(unique(df_metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # print(pal)

  # Run UMAP
  umap.out <- umap(df_PCA, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(df_PCA)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, df_metadata)

  # Plot UMAP
  plt <- ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(out_path,
                   paste(plot_name,meta_col,num_neighbor,"PCA.UMAP","png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)
}

make_PCA <- function(df_PCA, out_path,plot_name, meta_col){
  set.seed(123)

  # Make color palette
  n <- length(unique(df_metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # print(head(df_PCA))
  # print(tail(df_PCA))
  # print(dim(df_PCA))

  df_metadata <- as.data.frame(df_metadata) %>%
    filter(Run %in% rownames(df_PCA) )%>%
    tibble::column_to_rownames("Run")

  # print(head(df_metadata))
  # print(head(df_PCA))

  # Merge PCA results with metadata
  df_PCA <- data.frame(x = df_PCA[,1],  y = df_PCA[,2])
  pca.out.merge = cbind(df_PCA, df_metadata)

  # print(head(pca.out.merge))
  # print(tail(pca.out.merge))
  # print(dim(pca.out.merge))

  # Plot PCA
  plt <- ggplot(pca.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= "", sep = ' ')

  # Save plot
  ggsave(file.path(out_path,
                   paste(plot_name,meta_col,"PCA","png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)
}

df_to_UMAP <- function(input_df, output_dir, output_name){

  var <- apply(input_df[, -1], 1, var)
  param <- quantile(var, c(.1), na.rm=T)
  input_df_filt <- input_df[var > param & !is.na(var), ]

  # Transpose and format
  input_df_filt_t <- as.data.frame(t(input_df_filt))
  rownames(input_df_filt_t) <- colnames(input_df)

  print("Samples with NA, which will be dropped in PCA and UMAPs:")
  print(which(rowSums(is.na(input_df_filt_t))>0))

  # PCA
  prcomp.out = as.data.frame(prcomp(na.omit(input_df_filt_t), center=T,  scale = T)$x)
  print(dim(prcomp.out))
  
  # Plot PCAs
  for (meta in list("data_source", "LM22", "LM6", "sigil_general")){
    make_PCA(df_PCA = prcomp.out, out_path = paste0(output_dir, "/PCA/"), 
            plot_name = output_name, meta_col = paste0(meta))
    }


  # Plot variations of UMAPs with different numbers of neighbors
  lapply(c(20, 30), make_umap, meta_col="data_source",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)
  lapply(c(20, 30), make_umap, meta_col="LM22",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)
  lapply(c(20, 30), make_umap, meta_col="sigil_general",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)
  lapply(c(20, 30), make_umap, meta_col="LM6",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)

}

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
# scatter_per_gene <- function(df_long, ouput_prefix){

#   # Make scatter plot for each matching gene
#   ls_genes <- unique(df_long$gene)
#   for (g in ls_genes){
#     df_g <- df_long %>%
#       filter(gene == g)
#     p <- ggplot(aes(x=gene_z, y=splice_z, color = cell_ID, shape = event), data=df_g)+ 
#       geom_point() +
#       labs(title= g) + geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
#       geom_text(
#                 label= df_g$cell_types_splice,
#                 nudge_x = 0.05, nudge_y = 0.05,
#                 check_overlap =F, col = "black", size = 1
#               ) +
#       theme_minimal()

#     ggsave(plot = p, filename = paste0(ouput_prefix, g, ".png"))
#   }

# }
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#   #' Function to save pheatmaps to a pdf file
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }

# scatter_all <- function(df_long, output_path){
#   # Scatter plot with all matching events/genes
#   p <- ggplot(aes(x=gene_z, y=splice_z, color = gene), data=df_long)+ 
#     geom_point() +
#     geom_text(
#               label= df_long$gene,
#               nudge_x = 0.05, nudge_y = 0.05,
#               check_overlap =F, col = "black", size = 1.5
#             )
#   ggsave(plot = p, filename = output_path)
# }


upset_plot <- function(df_ref, val, output_name, ls_sets){
  # Build matrix
  #     cell1 cell2  cell3
  #event1 0      0     1
  #event2  1     1     0
  df_ref <- as.data.frame(df_ref)

  # Matrix filled of 0s
  df_upset<- data.frame(matrix(0,
                              ncol = length(unique(df_ref$cell_type)),
                              nrow = length(unique(df_ref[,paste0(val)]))
                              ))
  colnames(df_upset) <- unique(df_ref$cell_type)
  rownames(df_upset) <- unique(unique(df_ref[,paste0(val)]))

  # Loop through ref matrix rows and populate matrix
  for (row in 1:nrow(df_ref)) {
    cell <- as.character(df_ref[row, "cell_type"])
    event <- as.character(df_ref[row, paste0(val)])
    df_upset[event, cell] <- 1
  }

  if (ls_sets == "NA"){
    # nintersects 15 and nsets ncol(df_upset) takes 5 hrs to run 
      plotObject <- UpSetR::upset(df_upset, 
                              order.by = "freq",
                              keep.order = F, 
                              nintersects = 35, 
                              # nintersects = nrow(df_upset), 
                              # nsets= 2, 
                              nsets= ncol(df_upset), 
                              empty.intersections = "off")
    } else {
      plotObject <- UpSetR::upset(df_upset, 
                                order.by = "freq",
                                keep.order = F, 
                                nintersects = 15, 
                                sets = as.vector(unlist(ls_sets)), 
                                empty.intersections = "off")

    }
  
  pdf(file= paste0(output_name))
  print(plotObject)
  dev.off()


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
# if (!dir.exists(file.path(opt$out_dir,"/LM6_scatterplots/"))){
#   dir.create(file.path(opt$out_dir,"/LM6_scatterplots/"),
#               recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/dim_red/PCA/"))){
  dir.create(file.path(opt$out_dir,"/dim_red/PCA/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/dim_red/UMAP/"))){
  dir.create(file.path(opt$out_dir,"/dim_red/UMAP/"),
              recursive = TRUE, showWarnings = TRUE)}
# Read in files 

# Metadata
df_metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- df_metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Events in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "ref_matrix_PS/lm22_lm6_withinType_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 

# Merge cols to add LM tag (within type uses LM22 groups )
df_splice_ref <- df_splice_ref %>%
  mutate(cell_type_group = ifelse(group=="LM6",
                                paste(cell_type, "LM6", sep = "_"),
                                paste(cell_type, "LM22", sep = "_")))
print("splice ref mat")
print(head(df_splice_ref))
print(dim(df_splice_ref))
print(unique(df_splice_ref$group))

# Genes in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(
                              opt$geneDir,
                              "ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) %>%
              rename(gene = X)

# Merge cols to add LM tag (within type uses LM22 groups )
df_gene_ref <- df_gene_ref %>%
  mutate(cell_type_group = ifelse(group=="LM6",
                                paste(cell_type, "LM6", sep = "_"),
                                paste(cell_type, "LM22", sep = "_")))
print("gene ref mat")
print(head(df_gene_ref))
print(unique(df_gene_ref$group))

# Read in all gene exp 
df_exp <- read.table(file= paste0(
                              opt$geneDir,
                              "combined_kallisto_log2tpm_batch_corrected.csv"),
                    sep=",", header = TRUE, row.names=1)
df_exp <- df_exp %>% mutate_if(is.character,as.numeric)
print(head(df_exp))
print(dim(df_exp))

# Read in all MESA PS 
df_all_PS <- read.table(file = paste0(
                          opt$spliceDir,
                          "/batch_corr_mesa_allPS_LM22.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)

print(head(df_all_PS))
print(dim(df_all_PS))

# Read in MESA intron retention reference matrix
df_IR_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "ref_matrix_IR/lm22_lm6_withinType_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 
print(head(df_IR_ref))
print(dim(df_IR_ref))

# Merge cols to add LM tag (within type uses LM22 groups )
df_IR_ref <- df_IR_ref %>%
  mutate(cell_type_group = ifelse(group=="LM6",
                                paste(cell_type, "LM6", sep = "_"),
                                paste(cell_type, "LM22", sep = "_")))

# Read in MESA intron retention 
df_IR_table <- read.table(file = paste0(
                          opt$spliceDir,
                          "/batch_corr_mesa_ir_table_intron_retention_LM22.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_IR_table <- df_IR_table %>% mutate_if(is.character,as.numeric)

print(head(df_IR_table))
print(dim(df_IR_table))

######################
# PCA and UMAPS
#######################

# Gene UMAP
df_exp_ref <- df_exp %>%
  tibble::rownames_to_column('gene') %>%  
  filter(gene %in% df_gene_ref$gene) %>%
  tibble::column_to_rownames('gene')
# print(head(df_exp_ref))
print(dim(df_exp_ref))
df_to_UMAP(df_exp_ref, paste0(opt$out_dir, "/dim_red"), "gene_ref" )

# Splice UMAP
df_PS_ref <- df_all_PS %>%
  rownames_to_column("event") %>%
  filter(event %in% df_splice_ref$event) %>%
  column_to_rownames("event")
# print(head(df_PS_ref))
print(dim(df_PS_ref))
df_PS_ref_log <- as.data.frame(log2(df_PS_ref +1))
print(dim(df_PS_ref_log))
df_to_UMAP(df_PS_ref_log, paste0(opt$out_dir, "/dim_red"), "splice_ref" )

# IR UMAP
df_IR_table_ref <- df_IR_table %>%
  rownames_to_column("event") %>%
  filter(event %in% df_IR_ref$event) %>%
  column_to_rownames("event")
print(head(df_IR_table_ref))
print(dim(df_IR_table_ref))
df_IR_reflog <- as.data.frame(log2(df_IR_table_ref +1))
print(dim(df_IR_reflog))
df_to_UMAP(df_IR_reflog, paste0(opt$out_dir, "/dim_red"), "IR_ref" )

# Gene and splice UMAP
df_gene_ref_PS_ref <- rbind(df_exp_ref,df_PS_ref_log )
print(dim(df_PS_ref))
df_to_UMAP(df_gene_ref_PS_ref,paste0(opt$out_dir, "/dim_red"), "gene_and_splice_ref" )

# Gene and splice and IR UMAP
df_gene_ref_PS_ref_IR_ref <- rbind(df_gene_ref_PS_ref,df_IR_reflog )
print(dim(df_gene_ref_PS_ref_IR_ref))
df_to_UMAP(df_gene_ref_PS_ref_IR_ref,paste0(opt$out_dir, "/dim_red"), "gene_and_splice_and_IR_ref" )

##################################
# Formating Gene Exp Z-score df
#################################
# Gene exp : get medians and z-score across group
df_exp_LM22_med <- calcGroupMed(df_exp, unlist(unique(df_splice_ref$overlapping)), "LM22")
df_exp_LM22_med_z <- t(scale(t(df_exp_LM22_med)))
df_exp_LM6_med <- calcGroupMed(df_exp, unlist(unique(df_splice_ref$overlapping)), "LM6")
df_exp_LM6_med_z <- t(scale(t(df_exp_LM6_med)))

# Add LM6 / LM22 tags to column names because same cell name can be in both
colnames(df_exp_LM22_med_z) <- paste(colnames(df_exp_LM22_med_z),"LM22",sep=" ")
colnames(df_exp_LM6_med_z) <- paste(colnames(df_exp_LM6_med_z),"LM6",sep=" ")

# Gene exp: Combine LM6 and LM22 dfs
df_exp_LM6_LM22_med_z <- cbind(df_exp_LM6_med_z, df_exp_LM22_med_z)
print("Combined......")
print(head(df_exp_LM6_LM22_med_z))
print(dim(df_exp_LM6_LM22_med_z))

# Replace spaces in names with _ to match other df
str_conv_space <- function(in_str){
  paste(unlist(strsplit(as.character(in_str), split=" ")),
                                    collapse="_")
}

# Replace spaces with _ to match other df
new_exp_names <- unlist(lapply(colnames(df_exp_LM6_LM22_med_z),FUN= str_conv_space))
colnames(df_exp_LM6_LM22_med_z) <- new_exp_names
print(head(df_exp_LM6_LM22_med_z))

# Confirm new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type_group %in% colnames(df_exp_LM6_LM22_med_z))))
print("############################")

##################################
# Formating Splice Z-score df
#################################
df_lm22_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_PS/LM22_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1)


print(head(df_lm22_splice_z))
df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_PS/LM6_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1) 
print(head(df_lm6_splice_z))

# Add LM6 / LM22 tags to column names because same cell name can be in both
colnames(df_lm22_splice_z) <- paste(colnames(df_lm22_splice_z),"LM22",sep=" ")
colnames(df_lm6_splice_z) <- paste(colnames(df_lm6_splice_z),"LM6",sep=" ")

print(head(df_lm6_splice_z))
print(head(df_lm22_splice_z))

# Merge by row.names (cant cbind due to different number of events likely due to NANs)
df_splice_LM6_LM22_med_z <- merge(df_lm6_splice_z,df_lm22_splice_z,by=0) %>%
  column_to_rownames("Row.names")

# Splice: Combine LM6 and LM22
# df_splice_LM6_LM22_med_z <- cbind(df_lm6_splice_z, df_lm22_splice_z)
print(dim(df_lm6_splice_z))
print(dim(df_lm22_splice_z))
print(dim(df_splice_LM6_LM22_med_z))
print(head(df_splice_LM6_LM22_med_z))

# Replace spaces with _ to match other df
new_splice_names <- unlist(lapply(colnames(df_splice_LM6_LM22_med_z),FUN= str_conv_space))
colnames(df_splice_LM6_LM22_med_z) <- new_splice_names
print(head(df_splice_LM6_LM22_med_z))

# Check new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type_group %in% colnames(df_splice_LM6_LM22_med_z))))

##################################
# Formating IR Z-score df
#################################
df_lm22_IR_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_IR/LM22_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1)


print(head(df_lm22_IR_z))
df_lm6_IR_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_IR/LM6_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1) 
print(head(df_lm6_IR_z))

# Add LM6 / LM22 tags to column names because same cell name can be in both
colnames(df_lm22_IR_z) <- paste(colnames(df_lm22_IR_z),"LM22",sep=" ")
colnames(df_lm6_IR_z) <- paste(colnames(df_lm6_IR_z),"LM6",sep=" ")

print(head(df_lm6_IR_z))
print(head(df_lm22_IR_z))

# Merge by row.names (cant cbind due to different number of events likely due to NANs)
df_IR_LM6_LM22_med_z <- merge(df_lm6_IR_z,df_lm22_IR_z,by=0) %>%
  column_to_rownames("Row.names")

# IR: Combine LM6 and LM22
# df_IR_LM6_LM22_med_z <- cbind(df_lm6_IR_z, df_lm22_IR_z)
print(dim(df_lm6_IR_z))
print(dim(df_lm22_IR_z))
print(dim(df_IR_LM6_LM22_med_z))
print(head(df_IR_LM6_LM22_med_z))

# Replace spaces with _ to match other df
new_IR_names <- unlist(lapply(colnames(df_IR_LM6_LM22_med_z),FUN= str_conv_space))
colnames(df_IR_LM6_LM22_med_z) <- new_IR_names
print(head(df_IR_LM6_LM22_med_z))

# Check new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_IR_ref))
# print(head(df_IR_ref))
print(dim(df_IR_ref %>% filter(cell_type_group %in% colnames(df_IR_LM6_LM22_med_z))))

#################################################
# Adding all z-scores into df of splice ref events 
#################################################
# All genes that were in splice ref matrix
df_splice_ref_genes <- unique(df_splice_ref$overlapping)
print(length(df_splice_ref_genes))

# Add col to show if that gene was also in gene ref 
df_splice_ref <- df_splice_ref %>% 
  mutate(in_gene_ref = ifelse(overlapping %in% df_gene_ref$gene, paste(overlapping), ""))

df_splice_ref['gene_z'] <- df_splice_ref['splice_z'] <- df_splice_ref['IR_z'] <- NA

# Adding gene z_scores to splice ref df
for (gene in df_splice_ref_genes){
  #  For each gene get the cell_types the splice event was signifcant in
  df_  <- df_splice_ref %>% filter(overlapping == gene)

  # Check if gene in gene exp data (splicing can have lists for overlaps "AC129492.1,AC129492.4")
  if (gene %in% rownames(df_exp_LM6_LM22_med_z)){
    
    # Iterate over cell types and look up gene z score
    for (gene_cell in as.vector(unlist(unique(df_$cell_type_group)))){
      z <- df_exp_LM6_LM22_med_z[gene, gene_cell] 
      # Add gene z to each corresponding splice ref rows
      df_splice_ref <- df_splice_ref %>%
        mutate(gene_z = ifelse(((overlapping==gene) & (cell_type_group==gene_cell)), z, gene_z))
    }
  } else
  {print("no gene match")}
}

# Adding splice z_scores and IR z scores to splice ref df
for (this_event in df_splice_ref$event){
  df_  <- df_splice_ref %>% filter(event == this_event)
  # Iterate over cell types and look up Z scores 
  for (splice_cell in as.vector(unlist(unique(df_$cell_type_group)))){
    splice_z_val <- df_splice_LM6_LM22_med_z[this_event, splice_cell] 
    IR_z_val <- df_IR_LM6_LM22_med_z[this_event, splice_cell] 

    # Add zscores to ref df
    df_splice_ref <- df_splice_ref %>%
      mutate(splice_z = ifelse(((event==this_event) & (cell_type_group==splice_cell)), splice_z_val, splice_z)) %>%
      mutate(IR_z = ifelse(((event==this_event) & (cell_type_group==splice_cell)), IR_z_val, IR_z)) 
      
}
}

print(head(df_splice_ref))
print(unique(df_splice_ref$IR_z))

#################################################
# Adding all z-scores into df of IR ref events 
#################################################





################################################################################
# Quantify number of splice change with little gene change
###############################################################################
df_splice_ref$high_ratio <- NA
df_splice_ref <- df_splice_ref %>%
    mutate(ratio = abs(splice_z/gene_z)) %>%
    mutate(high_ratio = ifelse(ratio > 5,
                                paste(overlapping),
                                high_ratio)) %>%
    arrange(desc(ratio))

print("splice ref mat")
print(head(df_splice_ref, n = 50))
# print(tail(df_splice_ref, n = 50))

lowratio <- df_splice_ref %>%
  drop_na(ratio) %>%
  filter(ratio < 5) %>%
  arrange(ratio) %>%
  head(n = 100) 

# df_splice_ref %>%
#   filter(overlapping == "CSF2RA")

perc_event_high_ratio <- sum(!is.na(df_splice_ref$high_ratio))/nrow(df_splice_ref)
print("Percent of splice ref events with a splice z score to gene score ratio over 5")
print(perc_event_high_ratio)
for (i in unique(df_splice_ref$high_ratio)){
  cat(i)
  cat("\n")
}

for (i in unique(lowratio$overlapping)){
  cat(i)
  cat("\n")
}


#####################################################
# Scatter plot all splice events vs gene or IR
#####################################################
for (label_type in c("cell_type", "event", "overlapping", "group", "in_gene_ref", "high_ratio")){

  # Splice vs Gene________________________________________________________________
  p_vs_gene <- ggplot(aes(x=splice_z, y=gene_z, color = cell_type), data=df_splice_ref)+ 
    geom_point(size=.5) +
    labs(title= "") + 
    geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "grey") +
    geom_text(
              label= df_splice_ref[[label_type]],
              nudge_x = 0.05, nudge_y = 0.05,
              check_overlap =F, col = "black", size = 1
            ) +
    theme_minimal() +
    theme(legend.position="bottom", 
          legend.title = element_text(size = 5), 
          legend.text = element_text(size = 5),
          axis.title.x=element_text(size=7),
          axis.title.y=element_text(size=7))  +
    guides(color = guide_legend(nrow = 2))

  ggsave(plot = p_vs_gene, width = 12, height = 7, dpi = 400,
        filename = paste0(opt$out_dir, "/splice_ref_vs_gene_",label_type, ".png"))

  # Splice vs Gene Zoomed ______________________________________________________
  p_vs_gene_zoom <- ggplot(aes(x=splice_z, y=gene_z, color = cell_type), data=df_splice_ref)+ 
    geom_point(size=2) +
    labs(title= "") + 
    geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    geom_text(
              label= df_splice_ref[[label_type]],
              nudge_x = 0.01, nudge_y = 0.01,
              check_overlap =F, col = "black", size = 3
            ) +
    theme_minimal() +
    theme(legend.position="bottom", 
          legend.title = element_text(size = 5), 
          legend.text = element_text(size = 5),
          axis.title.x=element_text(size=7),
          axis.title.y=element_text(size=7))  +
    guides(color = guide_legend(nrow = 2)) +
    ylim(-1, 1) 

  ggsave(plot = p_vs_gene_zoom, width = 12, height = 7, dpi = 400,
        filename = paste0(opt$out_dir, "/splice_ref_vs_gene_",label_type, "_zoomed.png"))

  # Splice vs IR _____________________________________________________________
  p_vs_IR <- ggplot(aes(x=splice_z, y=IR_z, color = cell_type), data=df_splice_ref)+ 
    geom_point(size=.5) +
    labs(title= "") + 
    geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "grey") +
    geom_text(
              label= df_splice_ref[[label_type]],
              nudge_x = 0.05, nudge_y = 0.05,
              check_overlap =F, col = "black", size = 1
            ) +
    theme_minimal() +
    theme(legend.position="bottom", 
          legend.title = element_text(size = 5), 
          legend.text = element_text(size = 5),
          axis.title.x=element_text(size=7),
          axis.title.y=element_text(size=7))  +
    guides(color = guide_legend(nrow = 2))

  ggsave(plot = p_vs_IR, width = 12, height = 7, dpi = 400,
        filename = paste0(opt$out_dir, "/splice_ref_vs_IR_",label_type, ".png"))

}


# ########################################
# # Scatter for each cell type and group
# #######################################
# if (!dir.exists(file.path(opt$out_dir,"/splice_ref_vs_gene_by_cell/"))){
#   dir.create(file.path(opt$out_dir,"/splice_ref_vs_gene_by_cell/"),
#               recursive = TRUE, showWarnings = TRUE)}

# for (cell in unique(df_splice_ref$cell_type_group)){

#   df_splice_ref_subset <- df_splice_ref %>% 
#                           filter(cell_type_group == cell)

#   for (label_type in c( "event", "overlapping", "group", "in_gene_ref", "high_ratio")){

#     p <- ggplot(aes(x=splice_z, y=gene_z), data=df_splice_ref_subset)+ 
#       geom_point(size=2) +
#       labs(title= cell) + 
#       geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "black") +  
#       geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "black") +
#       geom_text(
#                 label= df_splice_ref_subset[[label_type]],
#                 nudge_x = 0.05, nudge_y = 0.05,
#                 check_overlap =F, col = "black", size = 2
#               ) +
#       theme_minimal() +
#       theme(legend.position="bottom", 
#             legend.title = element_text(size = 5), 
#             legend.text = element_text(size = 5),
#             axis.title.x = element_text(size=7),
#             axis.title.y = element_text(size=7))  +
#       guides(color = guide_legend(nrow = 2))

#     ggsave(plot = p, width = 12, height = 7, dpi = 400,
#           filename = paste0(opt$out_dir, "splice_ref_vs_gene_by_cell/splice_ref_vs_gene",cell,"_",label_type, ".png"))
#   }
# }


########################################
# Scatter plot all IR events 
#######################################











# #################
# # UpSet Plots
# ################
# T_cell_within_cell_types <-  c(
#     "T_cells_CD8",
#     "T_cells_CD4_naive",
#     #   "T cells CD4 memory resting",
#     #   "T cells CD4 memory  activated",
#     "T_cells_follicular_helper",
#     "T_cells_regulatory_(Tregs)",
#     "T_cells_gamma_delta")

# T_cell_within_cell_types <- c("T_cells_CD4_naive","T_cells_CD8", 
#     "T_cells_follicular_helper", "T_cells_gamma_delta"  )
     

# ls_upsets <- list(
#   # list(df_splice_ref, "event", paste0(opt$out_dir, "upsetplot_event2cell.pdf"), "NA"),
#   # list(df_splice_ref, "overlapping", paste0(opt$out_dir, "upsetplot_eventgene2cell.pdf"), "NA"),
#   # list(df_gene_ref, "gene", paste0(opt$out_dir, "upsetplot_DEG2cell.pdf", "NA")),

#   list(df_splice_ref, "event", paste0(opt$out_dir, "upsetplot_event2cell_Tcell.pdf"), T_cell_within_cell_types),
#   list(df_splice_ref, "overlapping", paste0(opt$out_dir, "upsetplot_eventgene2cell_Tcell.pdf"), T_cell_within_cell_types),
#   list(df_gene_ref, "gene", paste0(opt$out_dir, "upsetplot_DEG2cell_Tcell.pdf"), T_cell_within_cell_types)
# )

# foreach(i=ls_upsets, .packages=  c('magrittr', 'UpSetR')) %dopar% {
#   upset_plot(
#       df_ref = i[1],
#       val= i[2],
#       output_name = i[3], 
#       ls_sets = i[4]
#       )
#   }









##########################################
# OLD
##########################################
# # Map events to their signifcant cell type
# df_event2cell <- df_splice_ref %>%
#     group_by(event) %>%
#     summarize(context = list(cell_type)) %>%
#     mutate(cell_types = map_chr(context, toString)) %>%
#     select(event, cell_types) %>%
#     as.data.frame()
# print(head(df_event2cell))
# print(dim(df_event2cell))

# df_event2cell <- merge(x = df_event2cell, y = df_splice_ref, by="event") %>%
#   distinct(event, .keep_all = TRUE)

# print(head(df_event2cell))
# print(dim(df_event2cell))

# print("=====================================================================")

# # Genes in gene reference matrix 
# df_gene_ref <- read.csv(file= paste0(
#                               opt$geneDir,
#                               "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
#                          sep = "\t",header=TRUE) %>%
#               rename(gene = X)
# print("gene ref mat")
# print(head(df_gene_ref))

# print(unique(df_gene_ref$group))
# print(dim(df_gene_ref))


# # Map genes to to their signifcant cell type
# df_DEG2cell <- df_gene_ref %>%
#     group_by(gene) %>%
#     summarize(context = list(cell_type)) %>%
#     mutate(cell_types = map_chr(context, toString)) %>%
#     select(gene, cell_types) %>%
#     as.data.frame()

# print(head(df_DEG2cell))
# print(dim(df_DEG2cell))

# #######################
# # LM6 group z_scores
# ########################

# # Open LM6 splice z scores 
# df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,
#                                       "explore_ref_matrix/LM6_med_zscore.csv"),
#                             header=TRUE) %>%
#                     rename(event = X)

# # Add sig cell type and gene name
# df_lm6_splice_z <-  merge(x=df_lm6_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
#     select( -column_label2,-group, -column_label, -cell_type) %>%
#     rename(gene = overlapping)

# print("----------------------------------------")
# print(head(df_lm6_splice_z))
# print(dim(df_lm6_splice_z))

# # Open LM6 gene z scores 
# df_lm6_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE) %>%
#     rename(gene = X)
# # Add  cell_type to z scores
# df_lm6_gene_z <- merge(x=df_DEG2cell,y=df_lm6_gene_z,by="gene",all=TRUE) 
# print("df_lm6_gene_z")
# print(head(df_lm6_gene_z))
# print(dim(df_lm6_gene_z))


# # # Scatterplots 
# # df_lm6_merged_z_long <- format_merge(df_lm6_gene_z,df_lm6_splice_z )
# # print(head(df_lm6_merged_z_long))
# # scatter_all(df_lm6_merged_z_long, paste0(opt$out_dir, "/lm6_scatter_all.png"))
# # scatter_per_gene(df_lm6_merged_z_long, paste0(opt$out_dir, "/LM6_scatterplots/"))

# # #######################
# # # LM22
# # ########################
# # Open LM22 splice z scores 
# df_lm22_splice_z <- read.csv(file= paste0(opt$spliceDir,
#                                       "explore_ref_matrix/LM22_med_zscore.csv"),
#                             header=TRUE) %>%
#                     rename(event = X)

# # Add sig cell type and gene name
# df_lm22_splice_z <-  merge(x=df_lm22_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
#     select( -column_label2,-group, -column_label, -cell_type) %>%
#     rename(gene = overlapping)

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

