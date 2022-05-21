#!/usr/bin/env Rscript

library(optparse)
library(RColorBrewer)
library(foreach)
library(doParallel)

# Arguments
option_list <- list(

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
    c("-p", "--public_path"),
    type = "character",
    default = NULL,
    help = "public path to dir where the bigwigs will exist")

  )

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
# cl <- makeCluster(detectCores() - 2, outfile = "")
# registerDoParallel(cl)

colors <- c(
  #  14 Song
  # "139,0,0", #	dark red	#8B0000	(139,0,0)
  # "165,42,42", #brown	#A52A2A	(165,42,42)
 	# "178,34,34",#	firebrick	#B22222	(178,34,34)
 	# "220,20,60", #crimson	#DC143C	(220,20,60)
 	# "255,0,0", #red	#FF0000	(255,0,0)
 	# "255,99,71", #tomato	#FF6347	(255,99,71)
 	# "255,127,80", # coral	#FF7F50	(255,127,80)
 	# "205,92,92", #indian red	#CD5C5C	(205,92,92)
 	# "240,128,128", #light coral	#F08080	(240,128,128)
 	# "233,150,122", #dark salmon	#E9967A	(233,150,122)
  # "250,128,114",	#salmon	#FA8072	(250,128,114)
 	# "255,160,122", #light salmon	#FFA07A	(255,160,122)
 	# "255,69,0", #orange red	#FF4500	(255,69,0)
 	# "255,140,0") #dark orange	#FF8C00	(255,140,0))
 	


  #  orange	#FFA500	(255,165,0)


  # #  12 Choi
 	# "255,215,0", #gold	#FFD700	(255,215,0)
 	# "184,134,11", #dark golden rod	#B8860B	(184,134,11)
 	# "218,165,32", #golden rod	#DAA520	(218,165,32)
 	# "238,232,170", #pale golden rod	#EEE8AA	(238,232,170)
 	# "189,183,107", #dark khaki	#BDB76B	(189,183,107)
 	# "240,230,140", #khaki	#F0E68C	(240,230,140)
 	# "128,128,0", #olive	#808000	(128,128,0)
 	# "255,255,0", #yellow	#FFFF00	(255,255,0)
 	# "154,205,50", #yellow green	#9ACD32	(154,205,50)
 	# "85,107,47", #dark olive green	#556B2F	(85,107,47)
 	# "107,142,35", #olive drab	#6B8E23	(107,142,35)
 	# "124,252,0")
   
  
  #  # 50 Calderon 
  
  # "0,255,255",		#00FFFF		Cyan		87		180	100	50		
  # "155,221,255",		#9BDDFF		Columbia Blue		83		200	100	80		
  # "8,232,222",		#08E8DE		Bright Turquoise		79		177	93	47		
  # "135,206,235",		#87CEEB		Sky Blue		77		197	71	73		
  # "174,198,207",		#AEC6CF		Pastel Blue		76		196	26	75		
  # "48,213,200",		#30D5C8		Turquoise		73		175	66	51		
  # "0,123,167",		#007BA7		Cerulean		46		196	100	33		
  # "0,128,128",		#008080		Teal		43		180	100	25		
  # "47,79,79",		#2F4F4F		Dark Slate Gray		28		180	25	25		

  # "100,149,237",		#6495ED		Cornflower Blue		62		219	79	66
  # "0,127,255",		#007FFF		Azure		57		210	100	50
  # "3,138,168",		#5D8AA8		Air Force Blue		53		204	30	51
  # "0,71,171",		#0047AB		Cobalt		36		215	100	34

  # "181,126,220",		#B57EDC		Lavender		64		275	57	68		
  # "143,0,255",		#8F00FF		Violet		53		274	100	50		
  # "105,53,156",		#69359C		Purple Heart		38		270	49	41		

  # "209,159,232",		#D19FE8		Bright Ube		74		281	61	77		
  # "203,153,201",		#CB99C9		Pastel Violet		70		302	32	70		
  # "189,51,164",		#BD33A4		Byzantine		52		311	58	47		
  # "148,0,211",		#9400D3		Dark Violet		48		282	100	41		
  # "142,69,133",		#8E4585		Plum		43		307	35	41		
  # "139,0,139",		#8B008B		Dark Magenta		39		300	100	27		
  # "112,41,99",		#702963		Byzantium		32		311	46	30		

  # "251,204,231",		#FBCCE7		Classic Rose		88		326	85	89		
  # "244,154,194",		#F49AC2		Pastel Magenta		75		333	80	78		

  # "124,252,0",    #lawn green	#7CFC00	(124,252,0)
  # "127,255,0", #chartreuse	#7FFF00	(127,255,0)
 	
  # "173,255,47", # green yellow	#ADFF2F	(173,255,47)
 	# "0,100,0", # dark green	#006400	(0,100,0)
 	# "0,128,0", # green	#008000	(0,128,0)
 	# "34,139,34", # forest green	#228B22	(34,139,34)
 	# "0,255,0", # lime	#00FF00	(0,255,0)
 	# "50,205,50", # lime green	#32CD32	(50,205,50)
 	# "144,238,144", # light green	#90EE90	(144,238,144)
 	# "152,251,152", # pale green	#98FB98	(152,251,152)
 	# "143,188,143",  # dark sea green	#8FBC8F	(143,188,143)
 	# "0,250,154", # medium spring green	#00FA9A	(0,250,154)
 	# "0,255,127", # spring green	#00FF7F	(0,255,127)
 	# "46,139,87", # sea green	#2E8B57	(46,139,87)
 	# "102,205,170", # medium aqua marine	#66CDAA	(102,205,170)
 	# "60,179,113", # medium sea green	#3CB371	(60,179,113)
 	# "32,178,170", # light sea green	#20B2AA	(32,178,170)


  # "223,255,0",		#DFFF00		Chartreuse		90		68	100	50		
  # "191,255,0",		#BFFF00		Lime		87		75	100	50		
  # "167,252,0",		#A7FC00		Spring Bud		85		80	100	49		
  # "209,226,49",		#D1E231		Pear		81		66	75	54		
  # "141,182,0",		#8DB600		Apple Green		62		74	100	36		
  # "132,132,130",		#848482		Battleship Grey		52		60	1	51		
  # "128,128,0",		#808000		Olive		47		60	100	25		
  # "75,83,32"		#4B5320		Army Green		30		69	44	23		

  # Monaco 32 reusing colors
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
 	"255,140,0", #dark orange	#FF8C00	(255,140,0))
  
   	"255,215,0", #gold	#FFD700	(255,215,0)
 	"184,134,11", #dark golden rod	#B8860B	(184,134,11)
 	"218,165,32", #golden rod	#DAA520	(218,165,32)
 	"238,232,170", #pale golden rod	#EEE8AA	(238,232,170)
 	"189,183,107", #dark khaki	#BDB76B	(189,183,107)
 	"240,230,140", #khaki	#F0E68C	(240,230,140)
 	"128,128,0", #olive	#808000	(128,128,0)
 	"255,255,0", #yellow	#FFFF00	(255,255,0)
 	"154,205,50", #yellow green	#9ACD32	(154,205,50)
 	"85,107,47", #dark olive green	#556B2F	(85,107,47)
 	"107,142,35", #olive drab	#6B8E23	(107,142,35)
 	"124,252,0",

  "173,255,47", # green yellow	#ADFF2F	(173,255,47)
 	"0,100,0", # dark green	#006400	(0,100,0)
 	"34,139,34", # forest green	#228B22	(34,139,34)
 	"152,251,152", # pale green	#98FB98	(152,251,152)
 	"143,188,143",  # dark sea green	#8FBC8F	(143,188,143)
 	"0,250,154") # medium spring green	#00FA9A	(0,250,154)

# Make output directories
if (!dir.exists(file.path(opt$out_dir,"/bigwigs/"))){
  dir.create(file.path(opt$out_dir,"/bigwigs/"),
              recursive = TRUE, showWarnings = TRUE)}
if (!dir.exists(file.path(opt$out_dir,"/bams_ucsc/"))){
  dir.create(file.path(opt$out_dir,"/bams_ucsc/"),
              recursive = TRUE, showWarnings = TRUE)}

df_metadata <- read.csv(file = opt$manifest)
print(head(df_metadata))

print(df_metadata$source_name)
print(length(unique(df_metadata$source_name)))
print(length(colors))

# Map colors to each label
df_colormap <- data.frame(track_label=unique(df_metadata$source_name),
                             RGB=colors)

print(df_colormap)

#Check its existence
if (file.exists(paste0(opt$out_dir,"/trackDb.txt"))) {
  #Delete file if it exists
  file.remove(paste0(opt$out_dir,"/trackDb.txt"))
}

outfile = file(paste0(opt$out_dir,"/trackDb.txt"), open = 'a') # open in “a”ppend mode


for (i in df_metadata$Run){

  source_name <- as.character(subset(df_metadata, Run ==i)$source_name)
  # treatment <- as.character(subset(df_metadata, Run ==i)$treatment)
  # label <- as.character(subset(df_metadata, Run ==i)$sigil_cell_type_treatment)
  # color <- as.character(subset(df_colormap, track_label ==label)$RGB)
  color <- as.character(subset(df_colormap, track_label ==source_name)$RGB)

  cat("track",i, '\n', file = outfile, sep = ' ', append = TRUE)
  cat("bigDataUrl", paste0(opt$public_path,i,".bw"), '\n', file = outfile, sep = ' ', append = TRUE)
  cat("shortLabel",source_name, '\n', file = outfile, sep = ' ', append = TRUE)
  # cat("longLabel", i,source_name, treatment,'\n', file = outfile, sep = ' ', append = TRUE)
  cat("longLabel", i,source_name,'\n', file = outfile, sep = ' ', append = TRUE)
  cat("type bigWig",'\n', file = outfile, sep = ' ', append = TRUE)
  cat("color",color, '\n', file = outfile, sep = ' ', append = TRUE)
  cat("autoScale on", '\n', file = outfile, sep = ' ', append = TRUE)
  cat("maxHeightPixels 50", '\n', file = outfile, sep = ' ', append = TRUE)

  cat('\n',  file = outfile, sep = ' ', append = TRUE)

  # cat(i, '\n', file = 'filename', sep = '', append = TRUE)


}


