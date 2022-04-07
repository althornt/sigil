#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(uwot)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(UpSetR)

cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)


########################
# Functions
#########################
format_merge <- function(df_gene,  df_splice){
  # Making long dfs
  df_gene_z_long <- df_gene %>%
      tidyr::gather(., key="cell_ID", value = gene_z, -c(gene, cell_types)) %>%
      as.data.frame()

  df_splice_z_long <- df_splice %>%
      tidyr::gather(., key="cell_ID", value = splice_z, -c(gene, event, cell_types)) %>%
      as.data.frame()

  # Merging the long dfs 
  df_merged_z_long <-  df_splice_z_long %>%
          inner_join(df_gene_z_long, by = c("gene","cell_ID"), suffix=c("_splice","_gene")) %>% 
          tidyr::drop_na()

  df_merged_z_long$splice_z = as.numeric(as.character(df_merged_z_long$splice_z))
  df_merged_z_long$gene_z = as.numeric(as.character(df_merged_z_long$gene_z))

  return(df_merged_z_long)

}
scatter_per_gene <- function(df_long, ouput_prefix){

  # Make scatter plot for each matching gene
  ls_genes <- unique(df_long$gene)
  for (g in ls_genes){
    df_g <- df_long %>%
      filter(gene == g)
    p <- ggplot(aes(x=gene_z, y=splice_z, color = cell_ID, shape = event), data=df_g)+ 
      geom_point() +
      labs(title= g) + geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
      geom_text(
                label= df_g$cell_types_splice,
                nudge_x = 0.05, nudge_y = 0.05,
                check_overlap =F, col = "black", size = 1
              ) +
      theme_minimal()

    ggsave(plot = p, filename = paste0(ouput_prefix, g, ".png"))
  }

}
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#   #' Function to save pheatmaps to a pdf file
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }

scatter_all <- function(df_long, output_path){
  # Scatter plot with all matching events/genes
  p <- ggplot(aes(x=gene_z, y=splice_z, color = gene), data=df_long)+ 
    geom_point() +
    geom_text(
              label= df_long$gene,
              nudge_x = 0.05, nudge_y = 0.05,
              check_overlap =F, col = "black", size = 1.5
            )
  ggsave(plot = p, filename = output_path)
}

calcGroupMed <- function(df_exp, ls_gene, LM_type){
  str_LM <- as.character(LM_type)
  meta_row <- df_sample_annotations[c(paste0(LM_type)),]

  # Add LM type metadata row to expression
  df_exp_meta <- rbind(df_exp, meta_row )
  rownames(df_exp_meta)[length(rownames(df_exp_meta))] <-paste0(LM_type) # name new row
  print(ls_gene)

  ls_gene <- intersect(ls_gene, rownames(df_exp))

  print(ls_gene)

  # transpose df, summarize to median of each group
  df_t_med <- df_exp_meta %>%
    rownames_to_column("row_name") %>%  
    filter((row_name %in% ls_gene) | (row_name %in% c(paste0(LM_type)))) %>%  
    gather(var, value, -row_name) %>% 
    spread(row_name, value) %>%
    filter(get(LM_type) != "") %>%     #remove samples without a LM grouping
    select(c(paste0(str_LM),"var", ls_gene)) %>% #keep LM col and gene subset
    mutate_at(vars(as.vector(ls_gene)), as.numeric) %>% # convert to numeric
    group_by_at(LM_type) %>%
    summarise_at(vars(ls_gene), funs(median(., na.rm=TRUE))) %>%
    column_to_rownames(paste(str_LM)) %>%
    as.data.frame(.) #%>%
  #   # select_if(~ !any(is.na(.))) #remove NA

  df_med <-  df_t_med  %>%
    rownames_to_column("rowname") %>% 
    gather(var, value, -rowname)  %>% 
    spread(rowname, value) %>% 
    as.data.frame(.) %>%
    column_to_rownames("var")
    
  # #Remove genes with 0 variance which break scaling
  # df_med_var<- df_med[apply(df_med, 1, var) != 0, ]

  # heatmap_res <- pheatmap(
  #         main = paste0(" "),
  #         df_med_var,
  #         scale = "row",
  #         show_rownames=F,
  #         show_colnames=T,
  #         na_col = "grey"
  #         )

  # save_pheatmap_pdf(
  #         heatmap_res,
  #         paste0(opt$out_dir,"/",paste(LM_type),"_med_heatmap.pdf")
  #         )

  return(df_med)
}

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

print(opt$spliceDir)
print(opt$geneDir)

# Make output directories
# if (!dir.exists(file.path(opt$out_dir,"/LM6_scatterplots/"))){
#   dir.create(file.path(opt$out_dir,"/LM6_scatterplots/"),
#               recursive = TRUE, showWarnings = TRUE)}
# if (!dir.exists(file.path(opt$out_dir,"/LM22_scatterplots/"))){
#   dir.create(file.path(opt$out_dir,"/LM22_scatterplots/"),
#               recursive = TRUE, showWarnings = TRUE)}

# Read in files 

# Metadata
metadata <- read.csv(file = opt$metadata)
df_sample_annotations <- metadata %>%
  dplyr::select(Run,LM22,LM6, sigil_general, data_source) %>%
  tibble::column_to_rownames("Run") %>%
  t()

# Events in splice reference matrix 
df_splice_ref <- read.csv(file= paste0(
                                opt$spliceDir,
                                "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                          sep = "\t",header=TRUE) 
print("splice ref mat")
print(head(df_splice_ref))
print(dim(df_splice_ref))
print(unique(df_splice_ref$group))

# Genes in gene reference matrix 
df_gene_ref <- read.csv(file= paste0(
                              opt$geneDir,
                              "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
                         sep = "\t",header=TRUE) %>%
              rename(gene = X)
print("gene ref mat")
print(head(df_gene_ref))

#############
# MAIN PLOT: 
##############################
# LM22 splice ref vs all matched genes
##############################
# Adding Z-scores into the splice ref df 

# Merge cols to add LM tag (within type uses LM22 groups )
df_splice_ref <- df_splice_ref %>%
  mutate(cell_type_group = ifelse(group=="LM6",
                                paste(cell_type, "LM6", sep = "_"),
                                paste(cell_type, "LM22", sep = "_")))
print(head(df_splice_ref))

##################################
# Formating Gene Exp Z-score df
#################################

# Read in all gene exp 
df_exp <- read.table(file= paste0(
                              opt$geneDir,
                              "combined_kallisto_log2tpm_batch_corrected.csv"),
                    sep=",", header = TRUE, row.names=1)
df_exp <- df_exp %>% mutate_if(is.character,as.numeric)
print(head(df_exp))
print(dim(df_exp))

# Gene exp : get medians and z-score across group
df_exp_LM22_med <- calcGroupMed(df_exp, unlist(unique(df_splice_ref$overlapping)), "LM22")
df_exp_LM22_med_z <- t(scale(t(df_exp_LM22_med)))
df_exp_LM6_med <- calcGroupMed(df_exp, unlist(unique(df_splice_ref$overlapping)), "LM6")
df_exp_LM6_med_z <- t(scale(t(df_exp_LM6_med)))

# Add LM6 / LM22 tags to column names because same cell name can be in both
colnames(df_exp_LM22_med_z) <- paste(colnames(df_exp_LM22_med_z),"LM22",sep=" ")
colnames(df_exp_LM6_med_z) <- paste(colnames(df_exp_LM6_med_z),"LM6",sep=" ")

# Gene exp: Combine LM6 and LM22 dfs
df_exp_LM6_LM22_med_z <- cbind(df_exp_LM6_med_z, df_exp_LM22_med_z)
print("Combined......")
print(head(df_exp_LM6_LM22_med_z))
print(dim(df_exp_LM6_LM22_med_z))

# Replace spaces in names with _ to match other df
str_conv_space <- function(in_str){
  paste(unlist(strsplit(as.character(in_str), split=" ")),
                                    collapse="_")
}

# Replace spaces with _ to match other df
new_exp_names <- unlist(lapply(colnames(df_exp_LM6_LM22_med_z),FUN= str_conv_space))
colnames(df_exp_LM6_LM22_med_z) <- new_exp_names
print(head(df_exp_LM6_LM22_med_z))

# Confirm new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type_group %in% colnames(df_exp_LM6_LM22_med_z))))
print("############################")
##################################
# Formating Splice Z-score df
#################################
df_lm22_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix/LM22_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1)


print(head(df_lm22_splice_z))
df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,
                                      "explore_ref_matrix/LM6_med_zscore.csv"),
                            header=TRUE, check.names=FALSE, row.names = 1) 
print(head(df_lm6_splice_z))

# Add LM6 / LM22 tags to column names because same cell name can be in both
colnames(df_lm22_splice_z) <- paste(colnames(df_lm22_splice_z),"LM22",sep=" ")
colnames(df_lm6_splice_z) <- paste(colnames(df_lm6_splice_z),"LM6",sep=" ")

print(head(df_lm6_splice_z))
print(head(df_lm22_splice_z))


# Merge by row.names (cant cbind due to different number of events likely due to NANs)
df_splice_LM6_LM22_med_z <- merge(df_lm6_splice_z,df_lm22_splice_z,by=0) %>%
  column_to_rownames("Row.names")

# Splice: Combine LM6 and LM22
# df_splice_LM6_LM22_med_z <- cbind(df_lm6_splice_z, df_lm22_splice_z)
print(dim(df_lm6_splice_z))
print(dim(df_lm22_splice_z))
print(dim(df_splice_LM6_LM22_med_z))
print(head(df_splice_LM6_LM22_med_z))

# Replace spaces with _ to match other df
new_splice_names <- unlist(lapply(colnames(df_splice_LM6_LM22_med_z),FUN= str_conv_space))
colnames(df_splice_LM6_LM22_med_z) <- new_splice_names
print(head(df_splice_LM6_LM22_med_z))

# Check new columns str format matches ref mat rows
# Add stop if not equal 
print(dim(df_splice_ref))
print(dim(df_splice_ref %>% filter(cell_type_group %in% colnames(df_splice_LM6_LM22_med_z))))

#################################################
# Adding z-scores into df of splice ref events 
#################################################
print("###############################################################")
# All genes that were in splice ref matrix
df_splice_ref_genes <- unique(df_splice_ref$overlapping)
print(length(df_splice_ref_genes))

# Add col to show if that gene was also in gene ref 
df_splice_ref <- df_splice_ref %>% 
  mutate(in_gene_ref = ifelse(overlapping %in% df_gene_ref$gene, paste(overlapping), ""))

df_splice_ref['gene_z'] <- NA
df_splice_ref['splice_z'] <- NA
print(head(df_splice_ref))

print(head(df_exp_LM6_LM22_med_z))

# Adding gene z_scores to splice ref df
for (gene in df_splice_ref_genes){
  # print("start.............")
  # print(gene)
  #  For each gene get the cell_types the splice event was signifcant in
  df_  <- df_splice_ref %>% filter(overlapping == gene)
  # print(as.vector(unlist(unique(df_$cell_type))))
  # print(df_)

  # Check if gene in gene exp data (splicing can have lists for overlaps "AC129492.1,AC129492.4")
  if (gene %in% rownames(df_exp_LM6_LM22_med_z)){
    
    # Iterate over cell types
    for (gene_cell in as.vector(unlist(unique(df_$cell_type_group)))){
      # print(gene_cell)

      # Look up the gene exp val for gene and cell type
      z <- df_exp_LM6_LM22_med_z[gene, gene_cell] 
      # print(z)

      # Add gene z to each corresponding splice ref rows
      df_splice_ref <- df_splice_ref %>%
        mutate(gene_z = ifelse(((overlapping==gene) & (cell_type_group==gene_cell)), z, gene_z))
    }
  } else
  {print("no gene match")}
}


print(head(df_splice_ref))
# print(df_splice_ref)

# Adding splice z_scores to splice ref df
for (this_event in df_splice_ref$event){
  print("start.................")
  print(this_event)
  df_  <- df_splice_ref %>% filter(event == this_event)
  print(df_)

  # Iterate over cell types
  for (splice_cell in as.vector(unlist(unique(df_$cell_type_group)))){
    print(splice_cell)
    print(df_ %>% filter(cell_type_group ==splice_cell ))

    z <- df_splice_LM6_LM22_med_z[this_event, splice_cell] 
    print(z)

    print(df_splice_LM6_LM22_med_z %>% select(splice_cell) %>% filter(row.names(.)==this_event))

    df_splice_ref <- df_splice_ref %>%
       mutate(splice_z = ifelse(((event==this_event) & (cell_type_group==splice_cell)), z, splice_z))
}
}


print(head(df_splice_ref))

# ########################################
# # Scatter plot all
# #######################################

# for (label_type in c("cell_type", "event", "overlapping", "group", "in_gene_ref")){

#   p <- ggplot(aes(x=splice_z, y=gene_z, color = cell_type), data=df_splice_ref)+ 
#     geom_point(size=.5) +
#     labs(title= "") + 
#     geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "grey") +  
#     geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "grey") +
#     geom_text(
#               label= df_splice_ref[[label_type]],
#               nudge_x = 0.05, nudge_y = 0.05,
#               check_overlap =F, col = "black", size = 1
#             ) +
#     theme_minimal() +
#     theme(legend.position="bottom", 
#           legend.title = element_text(size = 5), 
#           legend.text = element_text(size = 5),
#           axis.title.x=element_text(size=7),
#           axis.title.y=element_text(size=7))  +
#     guides(color = guide_legend(nrow = 2))

#   ggsave(plot = p, width = 12, height = 7, dpi = 400,
#         filename = paste0(opt$out_dir, "/splice_ref_vs_gene_",label_type, ".png"))
# }


# ########################################
# # Scatter for each cell type and group
# #######################################
# if (!dir.exists(file.path(opt$out_dir,"/splice_ref_vs_gene_by_cell/"))){
#   dir.create(file.path(opt$out_dir,"/splice_ref_vs_gene_by_cell/"),
#               recursive = TRUE, showWarnings = TRUE)}

# for (cell in unique(df_splice_ref$cell_type_group)){

#   df_splice_ref_subset <- df_splice_ref %>% 
#                           filter(cell_type_group == cell)

#   for (label_type in c( "event", "overlapping", "group", "in_gene_ref")){

#     p <- ggplot(aes(x=splice_z, y=gene_z), data=df_splice_ref_subset)+ 
#       geom_point(size=2) +
#       labs(title= cell) + 
#       geom_hline(yintercept = 0, size = .25, linetype='dotted', color = "black") +  
#       geom_vline(xintercept = 0, size = .25, linetype='dotted', color = "black") +
#       geom_text(
#                 label= df_splice_ref_subset[[label_type]],
#                 nudge_x = 0.05, nudge_y = 0.05,
#                 check_overlap =F, col = "black", size = 2
#               ) +
#       theme_minimal() +
#       theme(legend.position="bottom", 
#             legend.title = element_text(size = 5), 
#             legend.text = element_text(size = 5),
#             axis.title.x = element_text(size=7),
#             axis.title.y = element_text(size=7))  +
#       guides(color = guide_legend(nrow = 2))

#     ggsave(plot = p, width = 12, height = 7, dpi = 400,
#           filename = paste0(opt$out_dir, "splice_ref_vs_gene_by_cell/splice_ref_vs_gene",cell,"_",label_type, ".png"))
#   }
# }

# ########################################
# # Zoom in on scatterplot of all
# #######################################
# for (label_type in c("cell_type", "event", "overlapping", "group", "in_gene_ref")){

#   p <- ggplot(aes(x=splice_z, y=gene_z, color = cell_type), data=df_splice_ref)+ 
#     geom_point(size=2) +
#     labs(title= "") + 
#     geom_hline(yintercept = 0, size = .5, linetype='dotted', color = "grey") +  
#     geom_vline(xintercept = 0, size = .5, linetype='dotted', color = "grey") +
#     geom_text(
#               label= df_splice_ref[[label_type]],
#               nudge_x = 0.01, nudge_y = 0.01,
#               check_overlap =F, col = "black", size = 3
#             ) +
#     theme_minimal() +
#     theme(legend.position="bottom", 
#           legend.title = element_text(size = 5), 
#           legend.text = element_text(size = 5),
#           axis.title.x=element_text(size=7),
#           axis.title.y=element_text(size=7))  +
#     guides(color = guide_legend(nrow = 2)) +
#     ylim(-.5, .5) 

#   ggsave(plot = p, width = 12, height = 7, dpi = 400,
#         filename = paste0(opt$out_dir, "/splice_ref_vs_gene_",label_type, "_zoomed.png"))
# }

##########################################
# UpSet plots - event vs cell type group
##########################################
# Build matrix
#     cell1 cell2  cell3
#event1 0      0     1
#event2  1     1     0
# print(head(df_splice_ref))

# # Matrix filled of 0s
# df_upset_event2cell <- data.frame(matrix(0,
#                             ncol = length(unique(df_splice_ref$cell_type_group)),
#                             nrow = length(unique(df_splice_ref$event))))
# colnames(df_upset_event2cell) <- unique(df_splice_ref$cell_type_group)
# rownames(df_upset_event2cell) <- unique(df_splice_ref$event)

# print(head(df_upset_event2cell))


# # Loop through ref matrix rows and populate matrix
# for (row in 1:nrow(df_splice_ref)) {
#   cell <- as.character(df_splice_ref[row, "cell_type_group"])
#   event <- as.character(df_splice_ref[row, "event"])
#   df_upset_event2cell[event, cell] <- 1
# }

# print(tail(df_upset_event2cell))


# plotObject <- UpSetR::upset(df_upset_event2cell, 
#                             order.by = "freq",
#                             keep.order = F, 
#                             nintersects = 15, 
#                             nsets= ncol(df_upset_event2cell), 
#                             empty.intersections = "off")

# pdf(file= paste0(opt$out_dir, "/upset_event2cell.pdf")) # or other device
# plotObject
# dev.off()


# print(dim(df_upset_event2cell))


# print(colSums(df_upset_event2cell))



upset_plot <- function(df_ref, val, output_name){
  # Build matrix
  #     cell1 cell2  cell3
  #event1 0      0     1
  #event2  1     1     0

  df_ref <- as.data.frame(df_ref)
  print(head(df_ref))
  print(val)

  # Matrix filled of 0s
  df_upset<- data.frame(matrix(0,
                              ncol = length(unique(df_ref$cell_type_group)),
                              nrow = length(unique(df_ref[,paste0(val)]))
                              ))
  colnames(df_upset) <- unique(df_ref$cell_type_group)
  rownames(df_upset) <- unique(unique(df_ref[,paste0(val)]))
  # print("---------------")
  print(head(df_upset))

  # Loop through ref matrix rows and populate matrix
  for (row in 1:nrow(df_ref)) {
    cell <- as.character(df_splice_ref[row, "cell_type_group"])
    event <- as.character(df_splice_ref[row, paste0(val)])
    df_upset[event, cell] <- 1
  }

  plotObject <- UpSetR::upset(df_upset, 
                            order.by = "freq",
                            keep.order = F, 
                            nintersects = 12, 
                            nsets= 10, 
                            nsets= ncol(df_upset), 
                            empty.intersections = "off")
  pdf(file= paste0(output_name))
  print(plotObject)
  dev.off()
  }


#################
# UpSet Plots
################
ls_upsets <- list(
  list(df_splice_ref, "event", paste0(opt$out_dir, "upset_plot_event2cell.pdf")),
  list(df_splice_ref, "overlapping", paste0(opt$out_dir, "upset_plot_eventgene2cell.pdf"))
  #list(df_gene_ref, gene, paste0(opt$out_dir, "/upset_DEG2cell.pdf"))
)

# Making each upset plot may be slow 
foreach(i=ls_upsets, .packages=  c('magrittr', 'UpSetR')) %dopar% {
  upset_plot(
      df_ref = i[1],
      val= i[2],
      output_name = i[3] 
      )
  }





#########################################################
# UpSet plots - Splice event's gene vs cell type group
# #########################################################

# # Matrix filled of 0s
# df_upset_gene2cell <- data.frame(matrix(0,
#                             ncol = length(unique(df_splice_ref$cell_type_group)),
#                             nrow = length(unique(df_splice_ref$overlapping))))
# colnames(df_upset_gene2cell) <- unique(df_splice_ref$cell_type_group)
# rownames(df_upset_gene2cell) <- unique(df_splice_ref$overlapping)

# print(head(df_upset_gene2cell))


# # Loop through ref matrix rows and populate matrix
# for (row in 1:nrow(df_splice_ref)) {
#   cell <- as.character(df_splice_ref[row, "cell_type_group"])
#   gene <- as.character(df_splice_ref[row, "overlapping"])
#   df_upset_gene2cell[gene, cell] <- 1
# }

# print(tail(df_upset_gene2cell))


# plotObject <- UpSetR::upset(df_upset_gene2cell, 
#                             order.by = "degree",
#                             keep.order = F, 
#                             nintersects = 15, 
#                             nsets= 17, 
#                             empty.intersections = "off")

# pdf(file= paste0(opt$out_dir, "/upset_gene2cell.pdf")) # or other device
# plotObject
# dev.off()

###########################################################
# UpSet plots - DEG vs cell type group
##########################################################







# # Map events to their signifcant cell type
# df_event2cell <- df_splice_ref %>%
#     group_by(event) %>%
#     summarize(context = list(cell_type)) %>%
#     mutate(cell_types = map_chr(context, toString)) %>%
#     select(event, cell_types) %>%
#     as.data.frame()
# print(head(df_event2cell))
# print(dim(df_event2cell))

# df_event2cell <- merge(x = df_event2cell, y = df_splice_ref, by="event") %>%
#   distinct(event, .keep_all = TRUE)

# print(head(df_event2cell))
# print(dim(df_event2cell))

# print("=====================================================================")

# # Genes in gene reference matrix 
# df_gene_ref <- read.csv(file= paste0(
#                               opt$geneDir,
#                               "/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv"),
#                          sep = "\t",header=TRUE) %>%
#               rename(gene = X)
# print("gene ref mat")
# print(head(df_gene_ref))

# print(unique(df_gene_ref$group))
# print(dim(df_gene_ref))


# # Map genes to to their signifcant cell type
# df_DEG2cell <- df_gene_ref %>%
#     group_by(gene) %>%
#     summarize(context = list(cell_type)) %>%
#     mutate(cell_types = map_chr(context, toString)) %>%
#     select(gene, cell_types) %>%
#     as.data.frame()

# print(head(df_DEG2cell))
# print(dim(df_DEG2cell))

# #######################
# # LM6 group z_scores
# ########################

# # Open LM6 splice z scores 
# df_lm6_splice_z <- read.csv(file= paste0(opt$spliceDir,
#                                       "explore_ref_matrix/LM6_med_zscore.csv"),
#                             header=TRUE) %>%
#                     rename(event = X)

# # Add sig cell type and gene name
# df_lm6_splice_z <-  merge(x=df_lm6_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
#     select( -column_label2,-group, -column_label, -cell_type) %>%
#     rename(gene = overlapping)

# print("----------------------------------------")
# print(head(df_lm6_splice_z))
# print(dim(df_lm6_splice_z))

# # Open LM6 gene z scores 
# df_lm6_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM6_med_zscore.csv"), header=TRUE) %>%
#     rename(gene = X)
# # Add  cell_type to z scores
# df_lm6_gene_z <- merge(x=df_DEG2cell,y=df_lm6_gene_z,by="gene",all=TRUE) 
# print("df_lm6_gene_z")
# print(head(df_lm6_gene_z))
# print(dim(df_lm6_gene_z))


# # # Scatterplots 
# # df_lm6_merged_z_long <- format_merge(df_lm6_gene_z,df_lm6_splice_z )
# # print(head(df_lm6_merged_z_long))
# # scatter_all(df_lm6_merged_z_long, paste0(opt$out_dir, "/lm6_scatter_all.png"))
# # scatter_per_gene(df_lm6_merged_z_long, paste0(opt$out_dir, "/LM6_scatterplots/"))

# # #######################
# # # LM22
# # ########################
# # Open LM22 splice z scores 
# df_lm22_splice_z <- read.csv(file= paste0(opt$spliceDir,
#                                       "explore_ref_matrix/LM22_med_zscore.csv"),
#                             header=TRUE) %>%
#                     rename(event = X)

# # Add sig cell type and gene name
# df_lm22_splice_z <-  merge(x=df_lm22_splice_z,y=df_event2cell,by="event",all=TRUE) %>%
#     select( -column_label2,-group, -column_label, -cell_type) %>%
#     rename(gene = overlapping)

# print("----------------------------------------")
# print(head(df_lm22_splice_z))
# print(dim(df_lm22_splice_z))

# # Open LM22 gene z scores 
# df_lm22_gene_z <- read.csv(file= paste0(opt$geneDir,"explore_ref_matrix/LM22_med_zscore.csv"), header=TRUE) %>%
#     rename(gene = X)
# # Add  cell_type to z scores
# df_lm22_gene_z <- merge(x=df_DEG2cell,y=df_lm22_gene_z,by="gene",all=TRUE) 
# print("df_lm22_gene_z")
# print(head(df_lm22_gene_z))
# print(dim(df_lm22_gene_z))

# # # Scatterplots 
# df_lm22_merged_z_long <- format_merge(df_lm22_gene_z,df_lm22_splice_z )
# scatter_all(df_lm22_merged_z_long, paste0(opt$out_dir, "/lm22_scatter_all.png"))
# scatter_per_gene(df_lm22_merged_z_long, paste0(opt$out_dir, "/LM22_scatterplots/"))

##############################
# Count splicing types 
##############################
	# - event - NA - no matched gene
	# - event - matched gene not in gene ref matrix
	# - event - matched gene in ref matrix same cell 
	# - event - matched gene in ref matrix different cell 
# print("Total number of events:")
# print(nrow(df_splice_ref))

# print("Total number of unique events:")
# print(length(unique(df_splice_ref$event)))

# print("Number of events with no gene:")
# print( nrow(df_splice_ref %>% 
#   filter(overlapping == "")))

# df_merged_z <- df_lm22_splice_z  %>%
#         inner_join(df_lm22_gene_z, by = c("gene"), suffix=c("_splice","_gene")) %>%
#         select(event, gene, cell_types_splice, gene,cell_types_gene ) %>%
#         arrange(gene)
        
# print("Number of events with matched gene in gene reference matrix:")
# print(nrow(df_merged_z))

# print("Number of spliced genes with matched gene in gene reference matrix:")
# print(length(unique(df_merged_z$gene)))

