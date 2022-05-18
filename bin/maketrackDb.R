#!/usr/bin/env Rscript

library(optparse)
library(RColorBrewer)
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
# cl <- makeCluster(detectCores() - 2, outfile = "")
# registerDoParallel(cl)

colors <- c(
  "139,0,0", #	dark red	#8B0000	(139,0,0)
  "165,42,42", #brown	#A52A2A	(165,42,42)
 	"178,34,34",#	firebrick	#B22222	(178,34,34)
 	"220,20,60", #crimson	#DC143C	(220,20,60)
 	"255,0,0", #red	#FF0000	(255,0,0)
 	"255,99,71", #tomato	#FF6347	(255,99,71)
 	"255,127,80", # coral	#FF7F50	(255,127,80)
 	"205,92,92", #indian red	#CD5C5C	(205,92,92)
 	"240,128,128", #light coral	#F08080	(240,128,128)
 	"233,150,122", #dark salmon	#E9967A	(233,150,122)
  "250,128,114",	#salmon	#FA8072	(250,128,114)
 	"255,160,122", #light salmon	#FFA07A	(255,160,122)
 	"255,69,0", #orange red	#FF4500	(255,69,0)
 	"255,140,0") #dark orange	#FF8C00	(255,140,0)
 	
  #  orange	#FFA500	(255,165,0)
 	# gold	#FFD700	(255,215,0)
 	# dark golden rod	#B8860B	(184,134,11)
 	# golden rod	#DAA520	(218,165,32)
 	# pale golden rod	#EEE8AA	(238,232,170)
 	# dark khaki	#BDB76B	(189,183,107)
 	# khaki	#F0E68C	(240,230,140)
 	# olive	#808000	(128,128,0)
 	# yellow	#FFFF00	(255,255,0)
 	# yellow green	#9ACD32	(154,205,50)
 	# dark olive green	#556B2F	(85,107,47)
 	# olive drab	#6B8E23	(107,142,35)
 	# lawn green	#7CFC00	(124,252,0)
 	# chartreuse	#7FFF00	(127,255,0)
 	# green yellow	#ADFF2F	(173,255,47)
 	# dark green	#006400	(0,100,0)
 	# green	#008000	(0,128,0)
 	# forest green	#228B22	(34,139,34)
 	# lime	#00FF00	(0,255,0)
 	# lime green	#32CD32	(50,205,50)
 	# light green	#90EE90	(144,238,144)
 	# pale green	#98FB98	(152,251,152)
 	# dark sea green	#8FBC8F	(143,188,143)
 	# medium spring green	#00FA9A	(0,250,154)
 	# spring green	#00FF7F	(0,255,127)
 	# sea green	#2E8B57	(46,139,87)
 	# medium aqua marine	#66CDAA	(102,205,170)
 	# medium sea green	#3CB371	(60,179,113)
 	# light sea green	#20B2AA	(32,178,170)


# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/bigwigs/"))){
  dir.create(file.path(opt$out_dir,"/bigwigs/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/bams_ucsc/"))){
  dir.create(file.path(opt$out_dir,"/bams_ucsc/"),
              recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = opt$manifest)
print(head(df_metadata))

print(length(unique(df_metadata$sigil_cell_type_treatment)))
print(df_metadata$sigil_cell_type_treatment)

print(colors)

# Map colors to each label
df_colormap <- data.frame(track_label=unique(df_metadata$sigil_cell_type_treatment),
                             RGB=colors)

print(df_colormap)

outfile = file(paste0(opt$out_dir,"/trackDb.txt"), open = 'a') # open in “a”ppend mode


public_path = "http://public.gi.ucsc.edu/~althornt/Song_tracks/"

for (i in df_metadata$Run){

  source_name <- as.character(subset(df_metadata, Run ==i)$source_name)
  treatment <- as.character(subset(df_metadata, Run ==i)$treatment)
  label <- as.character(subset(df_metadata, Run ==i)$sigil_cell_type_treatment)
  color <- as.character(subset(df_colormap, track_label ==label)$RGB)

  cat("track",i, '\n', file = outfile, sep = ' ', append = TRUE)
  cat("bigDataUrl", paste0(public_path,i,".bw"), '\n', file = outfile, sep = ' ', append = TRUE)
  cat("shortLabel",label, '\n', file = outfile, sep = ' ', append = TRUE)
  cat("longLabel", i,source_name, treatment,'\n', file = outfile, sep = ' ', append = TRUE)
  cat("type bigWig",'\n', file = outfile, sep = ' ', append = TRUE)
  cat("color",color, '\n', file = outfile, sep = ' ', append = TRUE)
  cat("autoScale on",color, '\n', file = outfile, sep = ' ', append = TRUE)

  cat('\n',  file = outfile, sep = ' ', append = TRUE)

  # cat(i, '\n', file = 'filename', sep = '', append = TRUE)



}



# foreach(i=df_metadata$Run) %dopar%{
#   # Find path to bam 
#   bam_path <- file.path(opt$bam_dir, paste0("STAR_", i),
#                   paste0(i,"Aligned.sortedByCoord.out.bam"))

#   # Remove extra chroms ; and update header 
#   cmd1 <- paste0("samtools view -b ", bam_path, 
#                 " 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT | samtools reheader ",
#                 opt$bam_header, " -  > ", opt$out_dir,"bams_ucsc/", i, ".bam")
#   print(cmd1)        
#   system(cmd1)        

#   # Make new index
#   cmd3 <- paste0("samtools index ",opt$out_dir,"bams_ucsc/", i, ".bam 2>&1")
#   print(cmd3)        
#   system(cmd3)        

#   # Make bigwigs from bam 
#   cmd4 <- paste0("bamCoverage -b ",opt$out_dir,"bams_ucsc/", i, ".bam -o ", 
#         opt$out_dir,"/bigwigs/",i,".bw --binSize 25 2>&1")
#   print(cmd4)       
#   system(cmd4)       
#   }
