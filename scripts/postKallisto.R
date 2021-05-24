
library(optparse)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(uwot)
library(ggplot2)
library(RColorBrewer)


set.seed(1234)

#arguments
option_list <- list(
  optparse::make_option(
    c("-i", "--kallisto_dir"),
    type = "character",
    default = NULL,
    help = "full path to directory with kallisto outputs"),

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
files <- file.path(opt$kallisto_dir, metadata$Run, "abundance.h5")

#transcripts to gene, used in tximport
edb <- EnsDb.Hsapiens.v86
tx2gene = transcripts(edb , columns=c("tx_id", "gene_name"),return.type="DataFrame")
all(file.exists(files))
names(files)<- metadata$Run
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene,
                         txOut = FALSE,ignoreTxVersion=TRUE,
                         countsFromAbundance = "scaledTPM")

#download counts
dat = txi.kallisto$counts
write.csv(dat,file.path(opt$out_dir,"combined_kallisto_tpm.csv"), row.names = TRUE)

#drop genes with low variance
getVar <- apply(dat[, -1], 1, var)
param <- median(getVar)
dat <- dat[getVar > param & !is.na(getVar), ]

#PCA
prcomp.out = as.data.frame(prcomp(as.data.frame(t(dat)), scale = F)$x)
prcomp.out$Run = rownames(prcomp.out)
prcomp.out.merge = merge(prcomp.out, y = metadata)

#UMAP function
make_umap <- function(num_neighbor,meta_col) {
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
          labs(title= paste("UMAP Kallisto: Cell types, n_neighbors =",num_neighbor, sep = ' '))

  #save
  ggsave(file.path(opt$out_dir,
          paste("kallisto_PCA_UMAP",meta_col,num_neighbor,"png", sep = '.')),
          device = "png",
          width = 12,
          dpi = 300)

}

#making variations of UMAPs
lapply(c(5,10,15,20,25,30), make_umap, meta_col="general")
lapply(c(5,10,15,20,25,30), make_umap, meta_col="Cell_type")
