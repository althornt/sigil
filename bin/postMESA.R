#!/usr/bin/env Rscript

library(optparse)
library(uwot)
library(ggplot2)
library(RColorBrewer)

#arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--mesa_PS"),
    type = "character",
    default = NULL,
    help = " `mesa_allPS.tsv` input file "),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"),

  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "full path to put outputs"))


#read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#make output directories
dir.create(file.path(opt$out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "UMAPs"), recursive = TRUE, showWarnings = FALSE)

#open files
metadata = read.csv(file = opt$metadata)
mesa_ps = read.csv(file = opt$mesa_PS, sep="\t", row.names = "cluster")

#keep junctions where less than 20% of samples are nans
cutoff = ncol(mesa_ps)*.20
mesa_ps_clean = mesa_ps[rowSums(is.na(mesa_ps))<cutoff,]

#replace remaining nan with row(junction) PS median
indx <- which(is.na(mesa_ps_clean), arr.ind = TRUE)
mesa_ps_clean[indx] <- apply(mesa_ps_clean, 1, median, na.rm = TRUE)[indx[,"row"]]

#drop genes with low variance
getVar <- apply(mesa_ps_clean[, -1], 1, var)
param <- median(getVar)
dat <- mesa_ps_clean[getVar > param & !is.na(getVar), ]


#UMAP function
make_umap <- function( num_neighbor,meta_col, df, tag) {

  set.seed(123)

  #PCA
  prcomp.out = as.data.frame(prcomp(as.data.frame(t(df)), scale = F)$x)
  prcomp.out$Run = rownames(prcomp.out)
  prcomp.out.merge = merge(prcomp.out, y = metadata)

  #make palette
  n <- length(unique(metadata[[meta_col]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  #umap
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  #merge
  umap.out.merge = merge(umap.out, metadata)

  #plot
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
          geom_point(size = 2) +
          theme_classic() +
          theme(legend.position="bottom",legend.title = element_blank()) +
          scale_color_manual(values=pal) +
          labs(title= paste("UMAP MESA: Cell types, n_neighbors =",num_neighbor,tag, sep = ' '))

  #save
  ggsave(file.path(opt$out_dir,
          paste("/UMAPs/MESA_PCA_UMAP",meta_col,num_neighbor,tag,"png", sep = '.')),
          device = "png",
          width = 12,
          dpi = 300)
}

#making variations of UMAPs
if("sigil_general" %in% colnames(metadata))
{
  lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="sigil_general", df = dat, tag="nan_adjusted")
  lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="sigil_general", df = na.omit(mesa_ps), tag="nan_dropped")
}

if("sigil_cell_type" %in% colnames(metadata))
{
  lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="sigil_cell_type", df = na.omit(mesa_ps), tag="nan_dropped")
  lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="sigil_cell_type", df = dat, tag="nan_adjusted")
}

if("sigil_cell_type_treatment" %in% colnames(metadata))
{
  lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="sigil_cell_type_treatment", df = na.omit(mesa_ps), tag="nan_dropped")
  lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="sigil_cell_type_treatment", df = dat, tag="nan_adjusted");
}


############################################
# NAN counting QC
############################################
#count nan per sample and per junction
column.nan.counts <- colSums(is.na(mesa_ps[,2:ncol(mesa_ps)]))
row.nan.counts <- rowSums(is.na(mesa_ps))
high_nan_samples = names(column.nan.counts[column.nan.counts > (dim(mesa_ps)[1]/2)])

#plot dist of nan per column
ggplot(NULL, aes(x=column.nan.counts)) + geom_histogram(binwidth = 1000)
ggsave(file.path(opt$out_dir,
        paste("qc/sample_nans", "png", sep = '.')),
        device = "png",
        width = 12,
        dpi = 300)

# plot dist of nan per junction
ggplot(NULL, aes(x=row.nan.counts)) + geom_histogram(binwidth = 2)
ggsave(file.path(opt$out_dir,
        paste("qc/junction_nans", "png", sep = '.')),
        device = "png",
        width = 12,
        dpi = 300)

#UMAP function labeling samples with high QC
make_umap_qc <- function( num_neighbor,meta_col, df, tag) {
  set.seed(123)

  #PCA
  prcomp.out = as.data.frame(prcomp(as.data.frame(t(df)), scale = F)$x)
  prcomp.out$Run = rownames(prcomp.out)
  prcomp.out.merge = merge(prcomp.out, y = metadata)

  #make palette
  n <- length(unique(metadata[[meta_col]]))
  # print(n)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  pal = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  #umap
  umap.out <- umap(prcomp.out, n_neighbors = num_neighbor, learning_rate = 0.5, init = "random")
  umap.out<- data.frame(x = umap.out[,1],  y = umap.out[,2])
  umap.out$Run <- rownames(prcomp.out)

  #merge
  umap.out.merge = merge(umap.out, metadata)

  #plot
  ggplot(umap.out.merge, aes(x, y, color = get(meta_col))) +
          geom_point(size = 2) +
          theme_classic() +
          theme(legend.position="bottom",legend.title = element_blank()) +
          scale_color_manual(values=pal) +
          labs(title= paste("UMAP MESA: Cell types, n_neighbors =",num_neighbor,tag, sep = ' '))+
          geom_text(aes(x, y, label = "NAN"), data = umap.out.merge[umap.out.merge$Run %in% high_nan_samples,],vjust = 0, nudge_y = 0.25)

  #save
  ggsave(file.path(opt$out_dir,
          paste("qc/MESA_PCA_UMAP",meta_col,num_neighbor,tag,"nanqc.png", sep = '.')),
          device = "png",
          width = 12,
          dpi = 300)

  }

#make UMAPS if there are samples with over half NAN
print("high_nan_samples:")
print(high_nan_samples)

if(length(high_nan_samples)>0)
{

  if("sigil_general" %in% colnames(metadata))
  {
    lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_general", df = na.omit(mesa_ps), tag="nan_dropped")
    lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_general", df = dat, tag="nan_adjusted")
  }

  if("sigil_cell_type" %in% colnames(metadata))
  {
    lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_cell_type", df = na.omit(mesa_ps), tag="nan_dropped")
    lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_cell_type", df = dat, tag="nan_adjusted")
  }

  if("sigil_cell_type_treatment" %in% colnames(metadata))
  {
    lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_cell_type_treatment", df = na.omit(mesa_ps), tag="nan_dropped")
    lapply(c(3,5,8,10,15,20,25), make_umap_qc, meta_col="sigil_cell_type_treatment", df = dat, tag="nan_adjusted")
  }
}
