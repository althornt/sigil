#!/usr/bin/env Rscript

library(optparse)
library(foreach)
library(doParallel)

# Arguments
option_list <- list(

  optparse::make_option(
    c("-d", "--bam_dir"),
    type = "character",
    default = NULL,
    help = "path to STAR directory containing STAR output directories per sample"),

optparse::make_option(
    c("-m", "--manifest"),
    type = "character",
    default = NULL,
    help = "path to manifest"),

  optparse::make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = NULL,
    help = "path to write output directory"),
  
  optparse::make_option(
    c("-b", "--bam_header"),
    type = "character",
    default = NULL,
    help = "new bam header")
  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
cl <- makeCluster(detectCores() - 2, outfile = "")
registerDoParallel(cl)


# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/bigwigs/"))){
  dir.create(file.path(opt$out_dir,"/bigwigs/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/bams_ucsc/"))){
  dir.create(file.path(opt$out_dir,"/bams_ucsc/"),
              recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = opt$manifest)
print(head(df_metadata))

foreach(i=df_metadata$Run) %dopar%{
  # Find path to bam 
  bam_path <- file.path(opt$bam_dir, paste0("STAR_", i),
                  paste0(i,"Aligned.sortedByCoord.out.bam"))

  # cmd0 <- paste0("samtools index ",bam_path)
  # print(cmd0)
  # system(cmd0)

  # Remove extra chroms ; and update header 
  cmd1 <- paste0("samtools view -b ", bam_path, 
                " 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT | samtools reheader ",
                opt$bam_header, " -  > ", opt$out_dir,"bams_ucsc/", i, ".bam")
  print(cmd1)        
  system(cmd1)        

  # Make new index
  cmd2 <- paste0("samtools index ",opt$out_dir,"bams_ucsc/", i, ".bam")
  print(cmd2)        
  system(cmd2)        

  # Make bigwigs from bam 
  cmd3 <- paste0("bamCoverage -b ",opt$out_dir,"bams_ucsc/", i, ".bam -o ", 
        opt$out_dir,"/bigwigs/",i,".bw --binSize 5 ")
  print(cmd3)       
  system(cmd3)       
  }
