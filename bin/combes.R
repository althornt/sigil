#!/usr/bin/env Rscript
library(optparse)
library(uwot)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(pheatmap)
library(purrr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)

###############
#  Functions
##############
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
  #' @param df_meta - df of metadata with Run, val, and data_source, main_label, group_label columns
  #' @param df_PS - df of MESA all PS file

  # Filter MESA all PS file to events of interest
  df_all_PS_sig_events <- df_PS %>%
    tibble::rownames_to_column('event') %>%
    dplyr::filter(event %in% ls_events) %>%
    dplyr::select(noquote(order(colnames(.)))) %>%
    tibble::column_to_rownames('event')

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

sample_cor <- function(df, tag){
    # Keep if junction number of nans is less than 20%
    cutoff = ncol(df)*.20
    df_clean = df[rowSums(is.na(df))<cutoff,]
    print(dim(df_clean))

    indx <- which(is.na(df_clean), arr.ind = TRUE)
    df_clean[indx] <- apply(df_clean, 1, median, na.rm = TRUE)[indx[,"row"]]

    print(dim(df_clean))


    cormat <- round(cor(df_clean, use = "everything"), 2)
    head(cormat)

    write.csv(cormat,paste0(opt$out_dir,tag,"_cormat.csv"), row.names = TRUE)

    # Get lower triangle of the correlation matrix
    get_lower_tri<-function(cormat){
        cormat[upper.tri(cormat)] <- NA
    return(cormat)
    }
    # Get upper triangle of the correlation matrix
    get_upper_tri <- function(cormat){
        cormat[lower.tri(cormat)]<- NA
    return(cormat)
    }

    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- melt(upper_tri, na.rm = TRUE)

    ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
        size = 5, hjust = 1))+
    theme(axis.text.y = element_text(size = 6)) +
    coord_fixed() 

    ggsave(paste0(opt$out_dir,tag,"_sample_cor.png" ))

}


######################
# Main 
######################

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

# Metadata
df_metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- df_metadata %>%
  dplyr::select(Run,main_label,group_label, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") 
#   %>%
#   t()

print(head(df_sample_annotations))
# Events in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "ref_matrix_PS/combined_main_group_withingroup_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 

# Read in all MESA PS 
df_all_PS <- read.table(file = paste0(
                          opt$spliceDir,
                          "/batch_corr_mesa_allPS.tsv"),
                          sep="\t", header = TRUE, row.names=1) 
df_all_PS <- df_all_PS %>% mutate_if(is.character,as.numeric)
df_PS_ref <- df_all_PS %>%
  tibble::rownames_to_column("junction") %>%
  filter(junction %in% df_splice_ref$event)

df_main_label_med <-  read.csv(file= paste0(
                                opt$spliceDir,
                                "explore_ref_matrix_PS/main_label_med.csv"),
                          sep = ",",header=TRUE, row.names = "X") 

# Read in Combes res
df_tcell_lung <- read.table("/mnt_/results/sigil_results_Combes_tcell_lung/mesa_out/mesa_allPS.tsv",  row.names = 1, header=T)
colnames(df_tcell_lung) <- paste("LungTcell", colnames(df_tcell_lung), sep = "_")

df_tcell_kidney <-  read.table("/mnt_/results/sigil_results_Combes_tcell_kidney/mesa_out/mesa_allPS.tsv", row.names = 1, header=T)
colnames(df_tcell_kidney) <- paste("KidneyTcell", colnames(df_tcell_kidney), sep = "_")

# Add Tcell lung to sample data
lungmeta <- data.frame(matrix(ncol = 4, nrow = ncol(df_tcell_lung)))
colnames(lungmeta) <- c("main_label", "group_label", "sigil_general", "data_source")
rownames(lungmeta) <- colnames(df_tcell_lung)
lungmeta[,"main_label"] <- "Lung Tcell"
lungmeta[,"group_label" ] <- "Lung Tcell"
lungmeta[,"sigil_general" ] <- "Lung Tcell"
lungmeta[,"data_source" ] <- "Combes"

print(head(lungmeta))

kidneymeta <- data.frame(matrix(ncol = 4, nrow = ncol(df_tcell_kidney)))
colnames(kidneymeta) <- c("main_label", "group_label", "sigil_general", "data_source")
rownames(kidneymeta) <- colnames(df_tcell_kidney)
kidneymeta[,"main_label"] <- "Kidney Tcell"
kidneymeta[,"group_label" ] <- "Kidney Tcell"
kidneymeta[,"sigil_general" ] <- "Kidney Tcell"
kidneymeta[,"data_source" ] <- "Combes"

# print(kidneymeta)
df_sample_annotations <- rbind(df_sample_annotations, lungmeta)
df_sample_annotations <- rbind(df_sample_annotations, kidneymeta)

ls_ref_juncs <- unique(df_splice_ref$event)
df_tcell_kidney_ref <- df_tcell_kidney %>%
    tibble::rownames_to_column("junction") %>%
    filter(junction %in% ls_ref_juncs)
print(dim(df_tcell_kidney_ref))

df_tcell_lung_ref <- df_tcell_lung %>%
    tibble::rownames_to_column("junction") %>%
    filter(junction %in% ls_ref_juncs)
    
print(dim(df_tcell_lung_ref))
print(head(df_tcell_lung_ref))
# merge tcell df 
df_tcell_kidney_lung_ref <- merge(df_tcell_lung_ref, df_tcell_kidney_ref, by="junction", all = T) %>%
    tibble::column_to_rownames("junction")
print(head(df_tcell_kidney_lung_ref))
print(dim(df_tcell_kidney_lung_ref))

write.csv(df_tcell_kidney_lung_ref,paste0(opt$out_dir,"combes_ref.csv"), row.names = TRUE)

# Get sample correlation using ref matrix junctions 
sample_cor(df_tcell_kidney_lung_ref, "combes_ref")


# Get sample correlation using random 500 junctions 
common <- intersect(rownames(df_tcell_lung), rownames(df_tcell_kidney))
sample <- sample(common, size = 500, replace = FALSE)
df_tcell_kidney_sample <- df_tcell_kidney %>%
    tibble::rownames_to_column("junction") %>%
    filter(junction %in% sample)

df_tcell_lung_sample <- df_tcell_lung %>%
    tibble::rownames_to_column("junction") %>%
    filter(junction %in% sample)
    
# merge tcell df 
df_tcell_kidney_lung_sample <- merge(df_tcell_lung_sample, df_tcell_kidney_sample, by="junction", all = T) %>%
    tibble::column_to_rownames("junction")

print(dim(df_tcell_kidney_lung_sample))
sample_cor(df_tcell_kidney_lung_sample, "combes_random500")


# Get sample correlation for ref matrix junctions in all samples 



df_tcell_kidney_lung_ref <- df_tcell_kidney_lung_ref %>%
      tibble::rownames_to_column("junction") 
print(head(df_PS_ref))
print(head(df_tcell_kidney_lung_ref))

df_tcell_kidney_lung_PS_ref<- merge(df_PS_ref, df_tcell_kidney_lung_ref, by="junction", all = T) %>%
    tibble::column_to_rownames("junction") %>%
    dplyr::select(noquote(order(colnames(.)))) 


print(dim(df_PS_ref))
print(dim(df_tcell_kidney_lung_ref))
print(dim(df_tcell_kidney_lung_PS_ref))

# sample_cor(df_tcell_kidney_lung_PS_ref, "healthy_combes_ref")

# Make heatmap with this cell types events only within this cell types samples
# make_pheatmap(df_tcell_kidney_lung_PS_ref$junction, paste0("healthy_combes"),
#         df_metadata, df_tcell_kidney_lung_PS_ref )

    #   stopifnot(rownames(df_sample_annotations) == colnames(df_all_PS_sig_events))


print(head(rownames(df_sample_annotations)))

print(head(colnames(df_tcell_kidney_lung_PS_ref)))

df_sample_annotations<- df_sample_annotations %>%
    select(group_label)

# make rowname 

# heatmap_res <- pheatmap(
#     main = paste0(" "),
#     df_tcell_kidney_lung_PS_ref,
#     # scale = "row",
#     show_rownames= F,
#     # show_colnames=T,
#     na_col = "grey",
#     annotation_col = df_sample_annotations,
#     clustering_method = "single")

# save_pheatmap_pdf(
# heatmap_res,
# paste0(opt$out_dir,"/heatmap.pdf"))




df_main_label_med <- df_main_label_med %>%
    tibble::rownames_to_column("junction") 
    # %>%
    # dplyr::select(noquote(order(colnames(.)))) 
# print(dim(df_tcell_kidney_lung_ref))
# print(head(df_tcell_kidney_lung_ref))
df_tcell_kidney_lung_ref <- df_tcell_kidney_lung_ref %>%
    tibble::column_to_rownames("junction")
    # %>%
    # as.numeric()

print(dim(df_tcell_kidney_lung_ref))
print(head(df_tcell_kidney_lung_ref))

# print(typeof(df_tcell_kidney_lung_ref))
# quit()

cutoff = ncol(df_tcell_kidney_lung_ref)*.30
df_tcell_kidney_lung_ref_clean = df_tcell_kidney_lung_ref[rowSums(is.na(df_tcell_kidney_lung_ref))<cutoff,]
print(dim(df_tcell_kidney_lung_ref_clean))

indx <- which(is.na(df_tcell_kidney_lung_ref_clean), arr.ind = TRUE)
df_tcell_kidney_lung_ref_clean[indx] <- apply(df_tcell_kidney_lung_ref_clean, 1, median, na.rm = TRUE)[indx[,"row"]]
print(warnings())
print(dim(df_tcell_kidney_lung_ref_clean))
# print(rowSums(is.na(df_tcell_kidney_lung_ref_drop)))
# quit()

print(head(df_main_label_med))
print(head(df_tcell_kidney_lung_ref_clean))
df_tcell_kidney_lung_ref_clean <- df_tcell_kidney_lung_ref_clean %>%
    tibble::rownames_to_column("junction")
df_main_label_med_tcell_kidney<- merge(df_main_label_med, df_tcell_kidney_lung_ref_clean, by="junction", all = F) %>%
    tibble::column_to_rownames("junction") 
    # %>%
    # as.numeric()
    
    #%>%
    #rop_na()
print(head(df_main_label_med_tcell_kidney))
print(tail(df_main_label_med_tcell_kidney))

print(dim(df_main_label_med_tcell_kidney))
# quit()

# print(tail(df_main_label_med_tcell_kidney))

# # print(dim(df_main_label_med_tcell_kidney))
# # quit()
# # heatmap_res_med <- pheatmap(
# #     main = paste0(" "),
# #     df_main_label_med_tcell_kidney,
# #     # scale = "row",
# #     show_rownames= F,
# #     show_colnames=T,
# #     na_col = "grey")

# # save_pheatmap_pdf(
# # heatmap_res_med,
# # paste0(opt$out_dir,"/heatmap_med.pdf"))

# scale rows
scaled_mat = t(scale(t(df_main_label_med_tcell_kidney)))

png(file=paste0(opt$out_dir,"/heatmap_med.png"),
    width = 5,
    height    = 3.25,
    units     = "in",
    res       = 1200)

ht <- ComplexHeatmap::Heatmap(scaled_mat,
                                show_row_names= FALSE,
                                column_names_gp = grid::gpar(fontsize = 3))
draw(ht)
dev.off()

# Get sample correlation using ref matrix junctions 
cormat <- round(cor(df_main_label_med_tcell_kidney, use = "everything"), 2)
head(cormat)

library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
col_fun(seq(-3, 3))

png(file=paste0(opt$out_dir,"/heatmap_med_corr.png"),
    width = 5,
    height    = 3.25,
    units     = "in",
    res       = 1200)

htcorr <- ComplexHeatmap::Heatmap(cormat,
                                row_names_gp = grid::gpar(fontsize = 3),
                                column_names_gp = grid::gpar(fontsize = 3))

draw(htcorr)
dev.off()