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
# Function to save pheatmaps to a file
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
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

  df_all_PS_sig_events_mat <- as.data.frame(lapply(df_all_PS_sig_events,
                                          function(x) as.numeric(as.character(x)))) 
                                          
  rownames(df_all_PS_sig_events_mat) <- rownames(df_all_PS_sig_events)

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

      # stopifnot(rownames(df_sample_annotations) == colnames(df_all_PS_sig_events))

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
plot_event <- function(sig_event, cell_type, out_dir, str_label){
  #' Make jitter plot for a given event in all samples
  #'
  #' @param sig_event - string for the significant event to ploit
  #' @param cell_type - string of the cell type the event was significant in
  #' @param out_dir - path to output directory
  #' @return NA

  df <- df_log2tpm_meta %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(gene %in% list(paste0(sig_event), paste0(str_label), "data_source"))%>%
    t() %>%
    as.data.frame()

  df_ <- df # copy df
  colnames(df_) <- c( "exp", paste0(str_label), "data_source") #add column names from first row

  df_ <- df_[-1,] %>% # drop first row 
          # dplyr::filter(PS != "NaN") %>% #drop samples with Nan
          dplyr::filter(get(str_label) != "") #drop samaples with out the cell type label

  df_$exp <- as.numeric(as.character(df_$exp))

  p <- ggplot( df_, aes(x = get(str_label), y = exp, color=data_source)) +
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
        str_label=meta_col, out_dir = plot_out_dir)

  # Get top negative delta events
  df_topN_DEG_up_reg <- df_res %>%
      dplyr::filter(padj <= .01 ) %>%
      dplyr::filter(log2FoldChange > 0 ) %>%
      dplyr::arrange(padj) %>%
      head(topN) %>%
      select(X)

  # Make plots for top positive events
  lapply(df_topN_DEG_up_reg$X,  plot_event, cell_type = LM22_type,
        str_label=meta_col, out_dir = plot_out_dir)

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

import_deseq2_within<- function(ls_cell_types, topN,  label, deseq2_dir, meta_col){
  #' Import results from MESA compare_sample_sets runs within a broader cell type
  #' (e.g. within T-cells). Find the top signficant events, make event level
  #' plots, make heatmaps of the events. This function calls import_deseq2(),
  #' unpack_import_css_res(), and make_pheatmap()
  #'
  #' @param ls_cell_types - list of  cell-types that were compared in
  #' compareWithinType.R script
  #' @param topN - integer; how many of the top splicing events to use
  #' @param label - string to use to represent cell type in output files
  #' @return ls_top_events - list containing 3 list - top positive events, top
  #' negative events, and top negative and positive


  # Get output files from compareWithinType script
  ls_css_file_names <- list.files(deseq2_dir,pattern = ".csv")
  print(ls_css_file_names)
  print(label)


  # For input cell type list , convert to filename
  ls_cell_types_file <- c()
  for (val in ls_cell_types){
    new_val <- paste0(gsub(" ","_", val), ".csv")
    ls_cell_types_file <- append(ls_cell_types_file, new_val)
  }

  # Intersect with the files that exist (Not all will have a mesa css output )
  ls_css_file_names_cell_type  <- intersect(ls_css_file_names, ls_cell_types_file)

  #Import files, find top signficant events, plot each event
  ls_res <- lapply(
                    ls_css_file_names_cell_type,
                    topN=topN,
                    import_deseq2,
                    meta_col =meta_col,
                    plot_out_dir =  paste0(opt$out_dir,"/ref_matrix/within_group/"),
                    deseq2_dir =  deseq2_dir,
                    comparison_label = "within_group")

  # If only 2 CSS files, they should have identical resuls
  # (comparing A vs B then B) , so only return the results of one
  if (length(ls_css_file_names_cell_type) <= 2) {
      #only keep first result
      df_res <- ls_res[[1]]
  } else {
      df_res <- dplyr::bind_rows(ls_res)
  }

  ls_top_events <- df_res$X

  # Filter metadata
  df_metadata_subset <- df_metadata %>%
    dplyr::filter(main_label %in% ls_cell_types) %>%
    droplevels(.) %>%
    dplyr::arrange(Run)

  # Filter all PS
  df_log2tpm_batch_corrrected_subset <- df_log2tpm_batch_corrrected %>%
    dplyr::select(as.vector(unlist(df_metadata_subset$Run)))

  # Make heatmap with this cell types events only within this cell types samples
  make_pheatmap(ls_top_events, paste0(label, "_gene_heatmap"),
          df_metadata_subset, df_log2tpm_batch_corrrected_subset )

  # Make heatmap with this cell types events and ALL samples
  make_pheatmap(ls_top_events, paste0(label,"_gene_heatmap_all_samples"),
          df_metadata, df_log2tpm_batch_corrrected )

  return(df_res)
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
if (!dir.exists(file.path(opt$out_dir,"/ref_matrix/"))){
  dir.create(file.path(opt$out_dir,"/ref_matrix/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/volcano_plots/"))){
  dir.create(file.path(opt$out_dir,"/volcano_plots/"),
              recursive = TRUE, showWarnings = TRUE)}


# Open files
df_metadata <- read.csv(file = paste0(opt$in_dir,"/metadata.csv"), header=TRUE, stringsAsFactors = FALSE)
df_sample_annotations <- df_metadata %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()
df_log2tpm_batch_corrrected <- read.csv(
    file = paste0(opt$in_dir,"/combined_kallisto_log2tpm_batch_corrected.csv"),
    header=TRUE, row.names = 1) 
# merge  metadata  and gene expression 
df_log2tpm_meta <- rbind(df_log2tpm_batch_corrrected, df_sample_annotations)

########################################
# Import main_label 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_main_label_de_file_names <- list.files(
                                  paste0(opt$out_dir,"/compare_main_label/deseq2_outputs/"),
                                  pattern = ".csv")

ls_cell_types_file <- list()
for (val in unique(list(df_metadata$main_label))){
    new_val <- paste0(gsub(" ","_", val), ".csv")
    ls_cell_types_file <- append(ls_cell_types_file, new_val)
  }

ls_de_file_names_cell_type  <- unlist(intersect(ls_main_label_de_file_names, ls_cell_types_file))

# Import, find signficant events and plot each one
ls_main_label_res <- foreach(i=ls_de_file_names_cell_type,
                      .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_deseq2(
      filename = i,
      topN = 20,
      meta_col="main_label",
      plot_out_dir = paste0(opt$out_dir,"/ref_matrix/main_label/"),
      deseq2_dir=paste0(opt$out_dir,"/compare_main_label/deseq2_outputs/"),
      comparison_label = "main_label"
      )
  }

df_main_label_res <- dplyr::bind_rows(ls_main_label_res)
print(head(df_main_label_res))
print(dim(df_main_label_res))

# Make heatmap using top events
make_pheatmap(df_main_label_res$X, "main_label_diff_gene_heatmap", df_metadata, df_log2tpm_batch_corrrected )


########################################
# Import group_label 1 vs all comparisons
#######################################

# Get all outputs from compare sample sets 1 vs all comparisons
ls_group_label_de_file_names <- list.files(
                                  paste0(opt$out_dir,"/compare_group_label/deseq2_outputs/"),
                                  pattern = ".csv")

ls_group_label_cell_types_file <- list()
for (val in unique(list(df_metadata$group_label))){
    new_val <- paste0(gsub(" ","_", val), ".csv")
    ls_group_label_cell_types_file <- append(ls_group_label_cell_types_file, new_val)
  }

ls_de_file_names_cell_type_group_label  <- unlist(intersect(ls_group_label_de_file_names, ls_group_label_cell_types_file))

# Import, find signficant events and plot each one
ls_group_label_res <- foreach(i=ls_de_file_names_cell_type_group_label,
                      .packages=c('magrittr','dplyr','ggplot2')) %dopar% {
    import_deseq2(
      filename = i,
      topN = 20,
      meta_col="group_label",
      plot_out_dir = paste0(opt$out_dir,"/ref_matrix/group_label/"),
      deseq2_dir=paste0(opt$out_dir,"/compare_group_label/deseq2_outputs/"),
      comparison_label = "group_label"
      )
  }

df_group_label_res <- dplyr::bind_rows(ls_group_label_res)
print(head(df_group_label_res))
print(dim(df_group_label_res))

# Make heatmap using top events
make_pheatmap(df_group_label_res$X, "group_label_diff_gene_heatmap", df_metadata, df_log2tpm_batch_corrrected )

########################################
# Import within type comparisons
# #######################################
# Match main cell type to group labels 
ls_group_cell_types <- unlist(unique(df_metadata[["group_label"]]))
ls_main_cell_types <- unlist(unique(df_metadata[["main_label"]]))

group2main <- list()
for (i in ls_group_cell_types){
  print(i)

  m <- df_metadata %>%
      filter(group_label == i) 
  if (length(unique(m$main_label)) > 1) {
      group2main<- append(group2main,list(list(i, unique(m$main_label)) ))
  }
}

print(group2main)
print("ppppppppp")
# Run in parallel
ls_within_res <- foreach(i=group2main,
                        .packages=c('magrittr','dplyr','ggplot2','pheatmap')) %dopar% {

      import_deseq2_within(
        ls_cell_types = i[[2]],
        topN = 20,
        label = i[[1]],
        deseq2_dir = paste0(opt$out_dir, "/compare_within_group/deseq2_outputs/"),
        meta_col = "main_label")

      }

df_within_res <- dplyr::bind_rows(ls_within_res, .id = "column_label")
print("df_within_res")
print(head(df_within_res))
print(dim(df_within_res))

# # Make heatmap using top events
make_pheatmap(df_within_res$X, "within_group_gene_heatmap", df_metadata, df_log2tpm_batch_corrrected)

########################################
# Combined reference matrix
######################################
df_combined_res <- dplyr::bind_rows(list("group_label"  = df_group_label_res,
                                        "main_label"= df_main_label_res,
                                        "within_group" = df_within_res),
                                        .id = "column_label2")
print("df_combined_res")
print(head(df_combined_res))
print(dim(df_combined_res))

write.table(df_combined_res,
            file=paste0(opt$out_dir,"/ref_matrix/combined_main_group_withingroup_combinedRefMat.tsv"),
            sep = "\t",row.names = FALSE, quote=F)

# Make heatmap using top events
make_pheatmap(df_combined_res$X, "combined_main_group_withingroup_gene_heatmap", 
              df_metadata, df_log2tpm_batch_corrrected)
