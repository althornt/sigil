#!/usr/bin/env Rscript

library(optparse)
library(uwot)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(stringr)
library(pheatmap)
library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

##################
# Functions
##################
# read_in_deseq2 <- function(file_path, topN) {
#   df_res <- read.csv(file = file_path, header = TRUE)

#   # Get cell type name from file path
#   str_cell_type <- tools::file_path_sans_ext(basename(file_path))

#   # Make volcano plot
#   volcano <- ggplot(data=df_res, aes(x=log2FoldChange, y=-log10(padj))) +
#     geom_point() + theme_minimal() +
#     geom_hline(yintercept=-log10(0.05), col="red") +
#     labs(title=paste(str_cell_type))
#   ggsave(
#     plot = volcano,
#     filename = paste0(opt$in_dir,"/volcano_plots/",str_cell_type,".png")
#     )

#   # Top N DEG by p-value regardless of log2FC
#   ls_topN_DEG_by_pval <- df_res %>%
#     dplyr::filter(padj < .05) %>%
#     dplyr::arrange(padj) %>%
#     head(topN) %>%
#     dplyr::pull(X)

#   # Top N DEG ranked by positive log2FC
#   ls_topN_DEG_up_reg <- df_res %>%
#     dplyr::filter(padj < .05) %>%
#     dplyr::arrange(desc(log2FoldChange)) %>%
#     head(topN) %>%
#     dplyr::pull(X)

#   # Top N DEG ranked by negative log2FC
#   ls_topN_DEG_down_reg <- df_res %>%
#     dplyr::filter(padj < .05) %>%
#     dplyr::arrange(desc(log2FoldChange)) %>%
#     tail(topN) %>%
#     dplyr::pull(X)

#   return(list("ls_topN_DEG_by_pval"= ls_topN_DEG_by_pval,
#               "ls_topN_DEG_up_reg" = ls_topN_DEG_up_reg,
#               "ls_topN_DEG_down_reg" = ls_topN_DEG_down_reg))

# }


# Function to save pheatmaps to a file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

list2heatmap <- function(ls_genes, title, output_name ){
  df_DEG_log2tpm_batch_corrrected <- df_log2tpm_batch_corrrected %>%
    dplyr::filter(X %in% ls_genes) %>%
    tibble::column_to_rownames(var  = "X")

  df_sample_annotations <- df_metadata %>%
  # dplyr::select(Run,sigil_general, LM22) %>%
    dplyr::select(Run, LM22, data_source) %>%
    tibble::column_to_rownames("Run")

  up_deg_heatmap <- pheatmap(
    main = paste0(title),
    df_DEG_log2tpm_batch_corrrected,
    scale = "row",
    show_rownames=F,
    show_colnames=F,
    annotation_col=df_sample_annotations)

  save_pheatmap_pdf(
    up_deg_heatmap,
    paste0(opt$out_dir,output_name,".pdf"))

}

make_umap <- function(num_neighbor,meta_col,df,out_path) {

  set.seed(123)
  df <- df %>%
    dplyr::select(-X)

  # Drop genes with low variance.
  getVar <- apply(df[, -1], 1, var)
  param <- median(getVar)
  log2trans_dat_filt <- df[getVar > param & !is.na(getVar), ]

  # Transpose and format
  log2trans_dat_filt_t <- as.data.frame(t(log2trans_dat_filt)) %>%
    dplyr::filter(!row.names(.) %in% c("gene"))
  rownames(log2trans_dat_filt_t) <- colnames(log2trans_dat_filt)

  # PCA.
  prcomp.out = as.data.frame(prcomp(log2trans_dat_filt_t, scale = F)$x)
  prcomp.out$Run = rownames(log2trans_dat_filt_t)
  # prcomp.out.merge = merge(prcomp.out, y = df)

  # Make color palette
  n <- length(unique(df_metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Run UMAP
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  # Merge UMAP results with metadata
  umap.out.merge = merge(umap.out, df_metadata)

  # Plot UMAP
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=pal) +
    labs(title= paste("UMAP Kallisto: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  # Save UMAP plot
  ggsave(file.path(opt$out_dir,
                   paste(out_path,meta_col,num_neighbor,"png", sep = '.')),
         device = "png",
         width = 12,
         dpi = 300)

}
plot_event <- function(sig_event, cell_type, out_dir, LM_type){
  #' Make jitter plot for a given event in all samples
  #'
  #' @param sig_event - string for the significant event to ploit
  #' @param cell_type - string of the cell type the event was significant in
  #' @param out_dir - path to output directory
  #' @return NA

  print(head(df_log2tpm_meta))
  print(sig_event)

  df <- df_log2tpm_meta %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(gene %in% list(paste0(sig_event), paste0(LM_type), "data_source"))%>%
    t() %>%
    as.data.frame()

  print(df)
  df_ <- df # copy df
  print(head(df_))

  colnames(df_) <- c( "exp", paste0(LM_type), "data_source") #add column names from first row

  df_ <- df_[-1,] %>% # drop first row 
          # dplyr::filter(PS != "NaN") %>% #drop samples with Nan
          dplyr::filter(get(LM_type) != "") #drop samaples with out the cell type label

  df_$exp <- as.numeric(as.character(df_$exp))

  print(head(df_))

  p <- ggplot( df_, aes(x = get(LM_type), y = exp, color=data_source)) +
      # geom_violin() +
      geom_jitter(position=position_jitter(0.15), alpha = 0.5, size = 2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(legend.position = "right") +
      scale_y_continuous() +
      labs(x = "cell type")

  ggsave(plot = p, filename = paste0(out_dir,cell_type,"/",sig_event,".png"))

}
import_deseq2 <- function(filename, topN, plot_out_dir, deseq2_dir, meta_col, comparison_label){
  #' Import results from MESA compare sample set script to get the top N
  #' significant events into lists
  #' @param filename -
  #' @param topN -
  #' @param plot_out_dir - path to output directory
  #' @param path tom esa compare sample set output directory
  #' @return  list of top N positive evnets, top negative events, top N
  #' negative and top N positive

  # Filename to string
  LM22_type <- stringr::str_trim(substr(filename, 1, nchar(filename)-4))
  print(LM22_type)

  # Make output directories
  if (!dir.exists(paste0(plot_out_dir,LM22_type))){
    dir.create(paste0(plot_out_dir,LM22_type),
     recursive = TRUE, showWarnings = TRUE)
  }

  # Open deseq2
  df_res <- read.csv(file = paste0(deseq2_dir,filename), header = TRUE)

  # Get top negative delta events
  df_topN_DEG_down_reg <- df_res %>%
      dplyr::filter(padj <= .01 ) %>%
      dplyr::filter(log2FoldChange < 0 ) %>%
      dplyr::arrange(padj) %>%
      head(topN) %>%
      select(X)

  # Make plots for top negative events
  lapply(df_topN_DEG_down_reg$X,  plot_event, cell_type = LM22_type,
        LM_type=meta_col, out_dir = plot_out_dir)

  # Get top negative delta events
  df_topN_DEG_up_reg <- df_res %>%
      dplyr::filter(padj <= .01 ) %>%
      dplyr::filter(log2FoldChange > 0 ) %>%
      dplyr::arrange(padj) %>%
      head(topN) %>%
      select(X)

  # Make plots for top positive events
  lapply(df_topN_DEG_up_reg$X,  plot_event, cell_type = LM22_type,
        LM_type=meta_col, out_dir = plot_out_dir)

  # Combine lists
  top_sig_neg_and_pos <- unlist(list(df_topN_DEG_down_reg$X,
                                    df_topN_DEG_up_reg$X))

  # Volcano plot labeling top events and coloring 
  df_res$point_color <- NA
  df_res <- df_res %>%  
    dplyr::arrange(padj) %>%
    dplyr::mutate(gglabel = case_when(((padj <.000000000001) | ((padj < .01  ) & (abs(log2FoldChange) > 2 )))~ X)) %>%
    dplyr::mutate(point_color = case_when(X %in% as.vector(top_sig_neg_and_pos)  ~ "Gene in the reference matrix",
                                        TRUE ~ "Gene not in the reference matrix"))

  volcano <- ggplot(data=df_res, aes(x=log2FoldChange, y=-log10(padj), color = point_color)) +
    geom_point(alpha= .7, size = 1) +
    theme_minimal() +
    geom_hline(yintercept=-log10(0.01), col="red") +
    labs(title=paste(LM22_type)) +
    geom_text(
          label= df_res$gglabel,
          nudge_x = 0.05, nudge_y = 0.05,
          check_overlap =F, col = "darkgreen", size = 2
        ) +
    scale_color_manual(name = "",
      values = c("Gene in the reference matrix" = "red",
                "Gene not in the reference matrix" = "black"),
      labels = c("Gene in the reference matrix",
                "Gene not in the reference matrix")) +
    theme(legend.position="bottom")
  
  ggsave(
    plot = volcano,
    filename = paste0(plot_out_dir,LM22_type,".png")
    )

  df_top_sig_neg_and_pos <- rbind(df_topN_DEG_down_reg, df_topN_DEG_up_reg)
  df_top_sig_neg_and_pos$group <- comparison_label
  df_top_sig_neg_and_pos$cell_type <- LM22_type

  return(df_top_sig_neg_and_pos)
}

###################
# MAIN
###################
# Arguments
option_list <- list(

  optparse::make_option(
    c("-i", "--in_dir"),
    type = "character",
    default = NULL,
    help = "path to inputs from combineGeneDE.R script"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "path to write outputs")
  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/refMatrix/"))){
  dir.create(file.path(opt$out_dir,"/refMatrix/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/volcano_plots/"))){
  dir.create(file.path(opt$out_dir,"/volcano_plots/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/UMAPs_DEG/"))){
  dir.create(file.path(opt$out_dir,"/UMAPs_DEG/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/heatmaps_DEG/"))){
  dir.create(file.path(opt$out_dir,"/heatmaps_DEG/"),
              recursive = TRUE, showWarnings = TRUE)}

# Open files
df_metadata <- read.csv(file = paste0(opt$in_dir,"/metadata.csv"), header=TRUE)
df_sample_annotations <- df_metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()
df_log2tpm_batch_corrrected <- read.csv(
    file = paste0(opt$in_dir,"/combined_kallisto_log2tpm_batch_corrected.csv"),
    header=TRUE, row.names = 1) 
# merge  metadata  and gene expression 
df_log2tpm_meta <- rbind(df_log2tpm_batch_corrrected, df_sample_annotations)

########################################
# Import LM22 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_lm22_css_file_names <- list.files(
                                  paste0(opt$out_dir,"/compare_LM22/deseq2_outputs/"),
                                  pattern = " .csv")[1] #keep space to avoid import of sample tables

print(ls_lm22_css_file_names)

# Import, find signficant events and plot each one
ls_lm22_res <- foreach(i=ls_lm22_css_file_names,
                      .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_deseq2(
      filename = i,
      topN = 20,
      meta_col="LM22",
      plot_out_dir = paste0(opt$out_dir,"/ref_matrix/LM22/"),
      deseq2_dir=paste0(opt$out_dir,"/compare_LM22/deseq2_outputs/"),
      comparison_label = "LM22"
      )
  }

# df_lm22_res <- dplyr::bind_rows(ls_lm22_res)
# print(head(df_lm22_res))
# print(dim(df_lm22_res))

# # Make heatmap using top events
# make_pheatmap(df_lm22_res$event, "LM22_diff_splicing_heatmap", metadata, all_PS )

# # Locate the deseq2 outputs
# ls_deseq_paths = list.files(
#   path = paste0(opt$in_dir,"/deseq2_outputs"),
#   pattern = "*" , full.names = TRUE)
# stopifnot(length(unique(df_metadata$LM22))==length(ls_deseq_paths))

# # Import deseq2 outputs, make volcano, and get list of top up DEG
# ls_topN_results <-lapply(ls_deseq_paths, read_in_deseq2, topN=20)

# # Split DEG lists from read_in_deseq2 into different lists
# ls_res_topN_DEG_by_pval <- unlist(lapply(c(1:length(ls_topN_results)),
#     function(x) ls_topN_results[[x]]$ls_topN_DEG_by_pval))
# ls_res_topN_upDEG <- unlist(lapply(c(1:length(ls_topN_results)),
#     function(x) ls_topN_results[[x]]$ls_topN_DEG_up_reg))
# ls_res_topN_downDEG <- unlist(lapply(c(1:length(ls_topN_results)),
#     function(x) ls_topN_results[[x]]$ls_topN_DEG_down_reg))




# ##################################
# # UMAPs and heatmaps using DEG
# ##################################
# # Using all genes for reference
# lapply(c(25,30), make_umap, meta_col="LM22",
#     df = df_log2tpm_batch_corrrected, out_path = "UMAPs_DEG/all_genes")

# # Filter df to top DEG and UMAP and heatmap
# df_log2tpm_batch_corrrected_topDEG <- df_log2tpm_batch_corrrected %>%
#     dplyr::filter(X %in% ls_res_topN_DEG_by_pval)
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22",
#     df = df_log2tpm_batch_corrrected_topDEG, out_path = "UMAPs_DEG/top_DEG_by_pval")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source",
#     df = df_log2tpm_batch_corrrected_topDEG, out_path = "UMAPs_DEG/top_DEG_by_pval")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general",
#     df = df_log2tpm_batch_corrrected_topDEG, out_path = "UMAPs_DEG/top_DEG_by_pval")
# list2heatmap(ls_res_topN_DEG_by_pval, "All top DEG","/heatmaps_DEG/top_DEG_by_pval" )

# # Filter df to top upregulated DEG and UMAP and heatmap
# df_log2tpm_batch_corrrected_top_up_DEG <- df_log2tpm_batch_corrrected %>%
#     dplyr::filter(X %in% ls_res_topN_upDEG)
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="LM22",
#     df = df_log2tpm_batch_corrrected_top_up_DEG, out_path = "UMAPs_DEG/top_up_DEG")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="data_source",
#     df = df_log2tpm_batch_corrrected_top_up_DEG, out_path = "UMAPs_DEG/top_up_DEG")
# lapply(c(5,10,15,20,25,30), make_umap, meta_col="sigil_general",
#     df = df_log2tpm_batch_corrrected_top_up_DEG, out_path = "UMAPs_DEG/top_up_DEG")
# list2heatmap(ls_res_topN_upDEG, "All top upregulated DEG", "/heatmaps_DEG/top_up_DEG")

# # Filter df to top downregulated DEG and UMAP and heatmap
# df_log2tpm_batch_corrrected_top_down_DEG <- df_log2tpm_batch_corrrected %>%
#     dplyr::filter(X %in% ls_res_topN_downDEG)
# lapply(c(5,10,15,20,25,30,35), make_umap, meta_col="LM22",
#     df = df_log2tpm_batch_corrrected_top_down_DEG, out_path = "UMAPs_DEG/top_down_DEG")
# lapply(c(5,10,15,20,25,30,35), make_umap, meta_col="data_source",
#     df = df_log2tpm_batch_corrrected_top_down_DEG, out_path = "UMAPs_DEG/top_down_DEG")
# lapply(c(5,10,15,20,25,30,35), make_umap, meta_col="sigil_general",
#     df = df_log2tpm_batch_corrrected_top_down_DEG, out_path = "UMAPs_DEG/top_down_DEG")
# list2heatmap(ls_res_topN_downDEG, "All top downregulated DEG", "/heatmaps_DEG/top_down_DEG")
