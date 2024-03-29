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

  print(head(umap.out.merge))

  # Plot UMAP
  plt <- ggplot(umap.out.merge, aes(x, y, color = get(meta_col), shape=data_source )) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("Cell types, n_neighbors =",num_neighbor, sep = ' ')) 
    # +
    # geom_text(
    #         label= umap.out.merge$Run,
    #         # vjust="inward",hjust="inward",
    #         nudge_x = 0.05, nudge_y = 0.05,
    #         check_overlap =F, col = "darkgreen", size = 2
    #       )

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
  for (meta in list("data_source", "main_label", "group_label", "sigil_general")){
    make_PCA(df_PCA = prcomp.out, out_path = paste0(output_dir, "/PCA/"), 
            plot_name = output_name, meta_col = paste0(meta))
    }


  # Plot variations of UMAPs with different numbers of neighbors
  lapply(c(20, 30), make_umap, meta_col="data_source",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)
  lapply(c(20, 30), make_umap, meta_col="main_label",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)
  lapply(c(20, 30), make_umap, meta_col="sigil_general",
    df_PCA = prcomp.out, out_path = paste0(output_dir, "/UMAP/"), plot_name = output_name)
  lapply(c(20, 30), make_umap, meta_col="group_label",
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
                              # nsets= 15, 
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
    
  print(sum(is.na(df_med)))
  print(colSums(is.na(df_med)))
  print(dim(df_med))

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

add_z_to_ref <- function(ref_df){
  # Function= Add gene, IR, and Splice z-score to an IR or Splice ref df
  
  # All genes that were in splice ref matrix
  df_splice_ref_genes <- unique(ref_df$overlapping)
  print(length(df_splice_ref_genes))

  # Add col to show if that gene was also in gene ref 
  ref_df <- ref_df %>% 
    mutate(in_gene_ref = ifelse(overlapping %in% df_gene_ref$gene, paste(overlapping), ""))

  ref_df['gene_z'] <- ref_df['splice_z'] <- ref_df['IR_z'] <- NA

  # Adding gene z_scores to splice ref df
  for (gene in df_splice_ref_genes){
    #  For each gene get the cell_types the splice event was signifcant in
    df_  <- ref_df %>% filter(overlapping == gene)

    # Check if gene in gene exp data (splicing can have lists for overlaps "AC129492.1,AC129492.4")
    if (gene %in% rownames(df_exp_group_label_main_label_med_z)){
      
      # Iterate over cell types and look up gene z score
      for (gene_cell in as.vector(unlist(unique(df_$cell_type_group)))){
        z <- df_exp_group_label_main_label_med_z[gene, gene_cell] 
        # Add gene z to each corresponding splice ref rows
        ref_df <- ref_df %>%
          mutate(gene_z = ifelse(((overlapping==gene) & (cell_type_group==gene_cell)), z, gene_z))
      }
    } else {
      print("missing gene:")
      print(gene)
      print(df_)

    }
  }

  # Adding splice z_scores and IR z scores to splice ref df
  for (this_event in ref_df$event){
    print(this_event)
    df_  <- ref_df %>% filter(event == this_event)
    # Iterate over cell types and look up Z scores 
    for (splice_cell in as.vector(unlist(unique(df_$cell_type_group)))){
      splice_z_val <- df_splice_group_label_main_label_med_z[this_event, splice_cell] 
      IR_z_val <- df_IR_group_label_main_label_med_z[this_event, splice_cell] 

      print(IR_z_val)

      # Add zscores to ref df
      ref_df <- ref_df %>%
        mutate(splice_z = ifelse(((event==this_event) & (cell_type_group==splice_cell)), splice_z_val, splice_z)) %>%
        mutate(IR_z = ifelse(((event==this_event) & (cell_type_group==splice_cell)), IR_z_val, IR_z)) 
        
  }
  }

  return(ref_df)
  
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
ls_out_paths <- list("/dim_red/PCA/","/dim_red/UMAP/",
                    "/upset_plots/","/upset_plots/",
                    "/zscore_scatter_plots/")
for (path in ls_out_paths){
  if (!dir.exists(paste0(opt$out_dir,path))){
    dir.create(paste0(opt$out_dir,path),
    recursive = TRUE, showWarnings = TRUE)
}
}


# Read in files 

# Metadata
df_metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- df_metadata %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Events in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "ref_matrix_PS/combined_main_group_withingroup_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 

# Merge cols to add LM tag (within type uses main_label groups )
df_splice_ref <- df_splice_ref %>%
  mutate(cell_type_group = ifelse(group=="group_label",
                                paste(cell_type, "group_label", sep = "_"),
                                paste(cell_type, "main_label", sep = "_")))
print("splice ref mat")
print(head(df_splice_ref))
print(dim(df_splice_ref))
print(unique(df_splice_ref$group))

# Genes in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(
                              opt$geneDir,
                              "ref_matrix/combined_main_group_withingroup_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) %>%
              rename(gene = X)

# Merge cols to add LM tag (within type uses main_label groups )
df_gene_ref <- df_gene_ref %>%
  mutate(cell_type_group = ifelse(group=="group_label",
                                paste(cell_type, "group_label", sep = "_"),
                                paste(cell_type, "main_label", sep = "_")))
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
                          "/batch_corr_mesa_allPS.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)

print(head(df_all_PS))
print(dim(df_all_PS))

# Read in MESA intron retention reference matrix
df_IR_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "ref_matrix_IR/combined_main_group_withingroup_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 
print(head(df_IR_ref))
print(dim(df_IR_ref))

# Merge cols to add LM tag (within type uses main_label groups )
df_IR_ref <- df_IR_ref %>%
  mutate(cell_type_group = ifelse(group=="group_label",
                                paste(cell_type, "group_label", sep = "_"),
                                paste(cell_type, "main_label", sep = "_")))

# Read in MESA intron retention 
df_IR_table <- read.table(file = paste0(
                          opt$spliceDir,
                          "/batch_corr_mesa_ir_table_intron_retention.tsv"),
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

ls_genes_from_splice_or_IR <- unlist(list(unlist(unique(df_splice_ref$overlapping)), unlist(unique(df_IR_ref$overlapping))))
print(head(ls_genes_from_splice_or_IR))
print(length(ls_genes_from_splice_or_IR))

df_exp_main_label_med <- calcGroupMed(df_exp, unlist(ls_genes_from_splice_or_IR), "main_label")
df_exp_main_label_med_z <- t(scale(t(df_exp_main_label_med)))
df_exp_group_label_med <- calcGroupMed(df_exp, unlist(ls_genes_from_splice_or_IR), "group_label")
df_exp_group_label_med_z <- t(scale(t(df_exp_group_label_med)))

# Add group_label / main_label tags to column names because same cell name can be in both
colnames(df_exp_main_label_med_z) <- paste(colnames(df_exp_main_label_med_z),"main_label",sep=" ")
colnames(df_exp_group_label_med_z) <- paste(colnames(df_exp_group_label_med_z),"group_label",sep=" ")

# Gene exp: Combine group_label and main_label dfs
df_exp_group_label_main_label_med_z <- cbind(df_exp_group_label_med_z, df_exp_main_label_med_z)
print("Combined......")
print(head(df_exp_group_label_main_label_med_z))
print(dim(df_exp_group_label_main_label_med_z))


# Replace spaces in names with _ to match other df
str_conv_space <- function(in_str){
  paste(unlist(strsplit(as.character(in_str), split=" ")),
                                    collapse="_")
}

# Replace spaces with _ to match other df
new_exp_names <- unlist(lapply(colnames(df_exp_group_label_main_label_med_z),FUN= str_conv_space))
colnames(df_exp_group_label_main_label_med_z) <- new_exp_names
print(head(df_exp_group_label_main_label_med_z))

# Confirm new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type_group %in% colnames(df_exp_group_label_main_label_med_z))))
print("############################")

##################################
# Formating Splice Z-score df
#################################
df_main_label_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_PS/main_label_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1)


print(head(df_main_label_splice_z))
df_group_label_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_PS/group_label_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1) 
print(head(df_group_label_splice_z))

# Add group_label / main_label tags to column names because same cell name can be in both
colnames(df_main_label_splice_z) <- paste(colnames(df_main_label_splice_z),"main_label",sep=" ")
colnames(df_group_label_splice_z) <- paste(colnames(df_group_label_splice_z),"group_label",sep=" ")

print(head(df_group_label_splice_z))
print(head(df_main_label_splice_z))

# Merge by row.names (cant cbind due to different number of events likely due to NANs)
df_splice_group_label_main_label_med_z <- merge(df_group_label_splice_z,df_main_label_splice_z,by=0) %>%
  column_to_rownames("Row.names")

# Splice: Combine group_label and main_label
# df_splice_group_label_main_label_med_z <- cbind(df_group_label_splice_z, df_main_label_splice_z)
print(dim(df_group_label_splice_z))
print(dim(df_main_label_splice_z))
print(dim(df_splice_group_label_main_label_med_z))
print(head(df_splice_group_label_main_label_med_z))

# Replace spaces with _ to match other df
new_splice_names <- unlist(lapply(colnames(df_splice_group_label_main_label_med_z),FUN= str_conv_space))
colnames(df_splice_group_label_main_label_med_z) <- new_splice_names
print(head(df_splice_group_label_main_label_med_z))

# Check new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type_group %in% colnames(df_splice_group_label_main_label_med_z))))

##################################
# Formating IR Z-score df
#################################
df_main_label_IR_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_IR/main_label_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1)


print(head(df_main_label_IR_z))
df_group_label_IR_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix_IR/group_label_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1) 
print(head(df_group_label_IR_z))

# Add group_label / main_label tags to column names because same cell name can be in both
colnames(df_main_label_IR_z) <- paste(colnames(df_main_label_IR_z),"main_label",sep=" ")
colnames(df_group_label_IR_z) <- paste(colnames(df_group_label_IR_z),"group_label",sep=" ")

print(head(df_group_label_IR_z))
print(head(df_main_label_IR_z))

# Merge by row.names (cant cbind due to different number of events likely due to NANs)
df_IR_group_label_main_label_med_z <- merge(df_group_label_IR_z,df_main_label_IR_z,by=0) %>%
  column_to_rownames("Row.names")

# IR: Combine group_label and main_label
# df_IR_group_label_main_label_med_z <- cbind(df_group_label_IR_z, df_main_label_IR_z)
print(dim(df_group_label_IR_z))
print(dim(df_main_label_IR_z))
print(dim(df_IR_group_label_main_label_med_z))
print(head(df_IR_group_label_main_label_med_z))

# Replace spaces with _ to match other df
new_IR_names <- unlist(lapply(colnames(df_IR_group_label_main_label_med_z),FUN= str_conv_space))
colnames(df_IR_group_label_main_label_med_z) <- new_IR_names
print(head(df_IR_group_label_main_label_med_z))

# Check new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_IR_ref))
# print(head(df_IR_ref))
print(dim(df_IR_ref %>% filter(cell_type_group %in% colnames(df_IR_group_label_main_label_med_z))))

###########################################################################
# Adding all z-scores into df of splice ref events and IR ref events 
###########################################################################
df_splice_ref <- add_z_to_ref(ref_df = df_splice_ref )
df_IR_ref <- add_z_to_ref(ref_df = df_IR_ref )

print(head(df_splice_ref))
print(dim(df_splice_ref))
print(colSums(is.na(df_splice_ref)))
print(length(unique(df_splice_ref$event)))

print("------------------------------------")

print(head(df_IR_ref))
print(dim(df_IR_ref))
print(colSums(is.na(df_IR_ref)))
print(length(unique(df_IR_ref$event)))


# ################################################################################
# # Quantify number of splice change with little gene change
# ###############################################################################
df_splice_ref$high_ratio <- NA
df_splice_ref <- df_splice_ref %>%
    mutate(ratio = abs(splice_z/gene_z)) %>%
    mutate(high_ratio = ifelse(ratio > 4,
                                paste(overlapping),
                                high_ratio)) %>%
    arrange(desc(ratio))

print("splice ref mat")
print(head(df_splice_ref, n = 50))
# print(tail(df_splice_ref, n = 50))

# Hundred lowest ratio genes 
lowratio <- df_splice_ref %>%
  drop_na(ratio) %>%
  filter(ratio < 4) %>%
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

print(head(df_splice_ref))
write.csv(df_splice_ref,paste0(opt$out_dir, "/splice_ref_zscores.csv"), row.names = FALSE)


df_IR_ref$high_ratio <- NA
df_IR_ref <- df_IR_ref %>%
    mutate(ratio = abs(IR_z/gene_z)) %>%
    mutate(high_ratio = ifelse(ratio > 4,
                                paste(overlapping),
                                high_ratio)) %>%
    arrange(desc(ratio))
write.csv(df_IR_ref,paste0(opt$out_dir, "/IR_ref_zscores.csv"), row.names = FALSE)


####################################################
# Scatter plot all splice events vs gene or IR
####################################################
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
        filename = paste0(opt$out_dir, "/zscore_scatter_plots/splice_ref_vs_gene_",label_type, ".png"))

  # Splice vs Gene Zoomed ______________________________________________________
  p_vs_gene_zoom <- ggplot(aes(x=splice_z, y=gene_z, color = cell_type), data=df_splice_ref)+ 
    geom_point(size=2) +
    labs(title= "", y = "Gene Expression Z-score", x = "Percent Spliced Z-score") +
    geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
    geom_text(
              label= df_splice_ref[[label_type]],
              nudge_x = 0.01, nudge_y = 0.01,
              check_overlap =F, col = "black", size = 3
            ) +
    theme_minimal() +
    theme(legend.position="bottom", 
          legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6),
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=8))  +
    guides(color = guide_legend(nrow = 3)) +
    ylim(-1, 1) 

  ggsave(plot = p_vs_gene_zoom, width = 12, height = 7, dpi = 400,
        filename = paste0(opt$out_dir, "/zscore_scatter_plots/splice_ref_vs_gene_",label_type, "_zoomed.png"))

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
        filename = paste0(opt$out_dir, "/zscore_scatter_plots/splice_ref_vs_IR_",label_type, ".png"))

}


# ########################################
# # Scatter for each cell type and group
# #######################################
if (!dir.exists(file.path(opt$out_dir,"/zscore_scatter_plots/splice_ref_vs_gene_by_cell/"))){
  dir.create(file.path(opt$out_dir,"/zscore_scatter_plots/splice_ref_vs_gene_by_cell/"),
              recursive = TRUE, showWarnings = TRUE)}

for (cell in unique(df_splice_ref$cell_type_group)){

  df_splice_ref_subset <- df_splice_ref %>% 
                          filter(cell_type_group == cell)

  for (label_type in c( "event", "overlapping", "group", "in_gene_ref", "high_ratio")){

    p <- ggplot(aes(x=splice_z, y=gene_z), data=df_splice_ref_subset)+ 
      geom_point(size=2) +
      labs(title= cell) + 
      geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "black") +  
      geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "black") +
      geom_text(
                label= df_splice_ref_subset[[label_type]],
                nudge_x = 0.05, nudge_y = 0.05,
                check_overlap =F, col = "black", size = 2
              ) +
      theme_minimal() +
      theme(legend.position="bottom", 
            legend.title = element_text(size = 5), 
            legend.text = element_text(size = 5),
            axis.title.x = element_text(size=7),
            axis.title.y = element_text(size=7))  +
      guides(color = guide_legend(nrow = 2))

    ggsave(plot = p, width = 12, height = 7, dpi = 400,
          filename = paste0(opt$out_dir, "zscore_scatter_plots/splice_ref_vs_gene_by_cell/splice_ref_vs_gene",cell,"_",label_type, ".png"))
  }
}


########################################
# Scatter plot all IR events 
#######################################

for (label_type in c("cell_type", "event", "overlapping", "group")){

  # IR vs Gene________________________________________________________________
  p_vs_gene <- ggplot(aes(x=IR_z, y=gene_z, color = cell_type), data=df_IR_ref)+ 
    geom_point(size=.5) +
    labs(title= "") + 
    geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "grey") +
    geom_text(
              label= df_IR_ref[[label_type]],
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
        filename = paste0(opt$out_dir, "/zscore_scatter_plots/IR_ref_vs_gene_",label_type, ".png"))

  # Splice vs IR _____________________________________________________________
  p_vs_IR <- ggplot(aes(x=IR_z, y=splice_z, color = cell_type), data=df_IR_ref)+ 
    geom_point(size=.5) +
    labs(title= "") + 
    geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "grey") +  
    geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "grey") +
    geom_text(
              label= df_IR_ref[[label_type]],
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
        filename = paste0(opt$out_dir, "/zscore_scatter_plots/IR_ref_vs_splice_",label_type, ".png"))

}



# # #################
# # # UpSet Plots
# # ################
# # T_cell_within_cell_types <-  c(
# #     "T_cells_CD8",
# #     "T_cells_CD4_naive",
# #     #   "T cells CD4 memory resting",
# #     #   "T cells CD4 memory  activated",
# #     "T_cells_follicular_helper",
# #     "T_cells_regulatory_(Tregs)",
# #     "T_cells_gamma_delta")

# # T_cell_within_cell_types <- c("T_cells_CD4_naive","T_cells_CD8", 
# #     "T_cells_follicular_helper", "T_cells_gamma_delta"  )
     

# # ls_upsets <- list(
# #   # list(df_splice_ref, "event", paste0(opt$out_dir, "upset_plots/upsetplot_event2cell.pdf"), "NA"),
# #   # list(df_splice_ref, "overlapping", paste0(opt$out_dir, "upset_plots/upsetplot_eventgene2cell.pdf"), "NA"),
# #   # list(df_gene_ref, "gene", paste0(opt$out_dir, "upset_plots/upsetplot_DEG2cell.pdf"), "NA"),
# #   list(df_IR_ref, "overlapping", paste0(opt$out_dir, "upset_plots/upsetplot_IRgene2cell.pdf"), "NA")

# #   # list(df_splice_ref, "event", paste0(opt$out_dir, "upset_plots/upsetplot_event2cell_Tcell.pdf"), T_cell_within_cell_types),
# #   # list(df_splice_ref, "overlapping", paste0(opt$out_dir, "upset_plots/upsetplot_eventgene2cell_Tcell.pdf"), T_cell_within_cell_types),
# #   # list(df_gene_ref, "gene", paste0(opt$out_dir, "upset_plots/upsetplot_DEG2cell_Tcell.pdf"), T_cell_within_cell_types)
# # )

# # foreach(i=ls_upsets, .packages=  c('magrittr', 'UpSetR')) %dopar% {
# #   upset_plot(
# #       df_ref = i[1],
# #       val= i[2],
# #       output_name = i[3], 
# #       ls_sets = i[4]
# #       )
# #   }


# Upset plot by gene 
# Genes affected by splice, IR, and gene
ls_gene_ref_genes <- unique(df_gene_ref$gene)
ls_splice_ref_genes <- unique(df_splice_ref$overlapping)
ls_IR_ref_genes <- unique(df_IR_ref$overlapping)

print("Number of gene ref genes:")
print(length(ls_gene_ref_genes))
print("Number of splice ref genes:")
print(length(ls_splice_ref_genes))
print("Number of IR ref genes:")
print(length(ls_IR_ref_genes))

# listInput <- list(gene = ls_gene_ref_genes,
#                  splice = ls_splice_ref_genes, 
#                  IR = ls_IR_ref_genes)

# plotObject <- upset(fromList(listInput), order.by = "freq")
# pdf(file= paste0(opt$out_dir, "/upset_plots/upsetplot_ref_genes_vs_splice_vs_IR.pdf"))
# print(plotObject)
# dev.off()

# # ##################################
# # # Compare gene lists
# # #################################

# Find genes in all 3
ls_ref_genes_splice_IR <- Reduce(intersect, 
                          list(ls_gene_ref_genes,
                              ls_splice_ref_genes,
                              ls_IR_ref_genes))
print("Intersection of all:")
print(ls_ref_genes_splice_IR)

# Find genes in splicing but nor ene 
ls_splice_ref_not_gene <- sort(unique(ls_splice_ref_genes[! ls_splice_ref_genes %in% ls_gene_ref_genes]))
print("Splice not gene:")
print(length(ls_splice_ref_not_gene))

for (i in ls_splice_ref_not_gene){
  cat(i)
  cat("\n")
}


# Find genes in IR but not gene 
ls_IR_ref_not_gene <- sort(unique(ls_IR_ref_genes[! ls_IR_ref_genes %in% ls_gene_ref_genes]))

print("IR not gene:")
print(length(ls_IR_ref_not_gene))

for (i in ls_IR_ref_not_gene){
  cat(i)
  cat("\n")
}

# Compare splice only to high ratio
print("High ratio and splice only genes ")
print(length(Reduce(intersect, 
                          list(unique(df_splice_ref$high_ratio),
                              ls_splice_ref_genes))))