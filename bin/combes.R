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

###############
#  Functions
##############
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
  tibble::column_to_rownames("Run") %>%
  t()

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
  tibble::rownames_to_column("event") %>%
  filter(event %in% df_splice_ref$event) %>%
  tibble::column_to_rownames("event")


# Read in Combes res
df_tcell_lung <- read.table("/mnt_/results/sigil_results_Combes_tcell_lung/mesa_out/mesa_allPS.tsv",  row.names = 1, header=T)
colnames(df_tcell_lung) <- paste("LungTcell", colnames(df_tcell_lung), sep = "_")

df_tcell_kidney <-  read.table("/mnt_/results/sigil_results_Combes_tcell_kidney/mesa_out/mesa_allPS.tsv", row.names = 1, header=T)
colnames(df_tcell_kidney) <- paste("KidneyTcell", colnames(df_tcell_kidney), sep = "_")

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

