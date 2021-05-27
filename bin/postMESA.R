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

#make output directory
dir.create(opt$out_dir)

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
  print(n)
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
          paste("MESA_PCA_UMAP",meta_col,num_neighbor,tag,"png", sep = '.')),
          device = "png",
          width = 12,
          dpi = 300)

}

#making variations of UMAPs
lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="general", df = dat, tag="nan_adjusted")
lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="cell_type", df = dat, tag="nan_adjusted")
lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="general", df = na.omit(mesa_ps), tag="nan_dropped")
lapply(c(3,5,8,10,15,20,25), make_umap, meta_col="cell_type", df = na.omit(mesa_ps), tag="nan_dropped")
