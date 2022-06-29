#!/usr/bin/env Rscript
library(optparse)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(cowplot)

cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

sigil_out_path = "/mnt/results_sigil_combine/sigil_results_SongChoi_newlabels_20220614/"
fig_output_path = "/mnt/compute_overlaps/"
if (!dir.exists(fig_output_path)){
    dir.create(fig_output_path,
    recursive = TRUE, showWarnings = TRUE)
    }

# Read in sigil sets
df_splice_set <- read.csv(file = paste0(sigil_out_path,
                          "combine_mesa_out/splice_set_PS/df_splice_sets.csv"))
print(head(df_splice_set))
df_IR_set <- read.csv(file = paste0(sigil_out_path,
                          "combine_mesa_out/splice_set_IR/df_splice_sets.csv"))
print(head(df_IR_set))
df_gene_set <- read.csv(file = paste0(sigil_out_path,
                          "combine_gene_out/gene_sets/df_gene_sets.csv"))
print(head(df_gene_set))

df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"), stringsAsFactors=FALSE)
ls_group_cell_types <- unlist(unique(df_metadata[["group_label"]]))
ls_main_cell_types <- unlist(unique(df_metadata[["main_label"]]))


# print(head(comb))


# for (i in comb){
#     print(i)
# }


compare_sets <-  function(df_sets, tag,total,col){
    ls_phyper <- list()
    ls_sets <- unique(df_sets$set)

    for (setA in ls_sets) {
        for (setB in ls_sets) {
            #skip if the same set
            if (setA==setB) break

            if (setA %in% c("T_Cells_CD8:+_within_UP","T_Cells_CD4:+_within_DN" ) &
                    setB %in% c("T_Cells_CD8:+_within_UP","T_Cells_CD4:+_within_DN" )
            ) break

            if (setA %in% c("B_cells_naive_within_UP","B_cells_memory_within_DN" ) &
                    setB %in% c("B_cells_naive_within_UP","B_cells_memory_within_DN" )
            ) break

        #        [1] "T_Cells_CD8:+_within_UP : T_Cells_CD4:+_within_DN"
        # [1] "T_Cells_CD8:+_within_DN : T_Cells_CD4:+_within_UP"
        # [1] "B_cells_naive_within_UP : B_cells_memory_within_DN"
        # [1] "B_cells_naive_within_DN : B_cells_memory_within_UP

            ls_events_A <- df_sets %>% 
                filter(set==setA) %>%
                pull(col)
            
            #    print(length(ls_events_A))
            ls_events_B <- df_sets %>% 
                filter(set==setB) %>%
                pull(col)
            #    print(length(ls_events_B))

            overlap <- length(intersect(ls_events_A,ls_events_B))
                # print(overlap)

            #Run hypergeometric test to find enrichment
            pval <- phyper(
                            overlap-1, 
                            length(ls_events_B),
                            total-length(ls_events_B), 
                            length(ls_events_A),
                            lower.tail= FALSE
                            )
                # print(pval)
            
                ls_phyper[[paste(setA,setB,sep=" : ")]] <- pval

                # print(paste(setA,setB,sep=" : "))
            
                if (overlap ==200){print(paste(setA,setB,sep=" : "))}
            }
        }

        # Make df of pvalues
        df_phyper <- do.call("rbind",ls_phyper) %>% as.data.frame()
        colnames(df_phyper) <- c("pval")

        # Make histogram of pvalues
        hist <- ggplot(df_phyper, aes(x=pval)) + 
            geom_histogram() +
            ggtitle(tag) +
            ylim(0, 4000) +
            theme_classic() 
        ggsave(paste0(fig_output_path,"phyper_hist_",tag,".png"), width=10, height=5, unit="cm")

        return(list("res" = df_phyper, "hist"=hist))
}

# Run function to get overlap among sigil sets 
phyper_splice <- compare_sets(df_splice_set, "splice", length(unique(df_splice_set$event)), "event")
print(head(phyper_splice$res))
print(dim(phyper_splice$res))

phyper_IR <- compare_sets(df_IR_set, "IR",length(unique(df_IR_set$event)), "event")
print(head(phyper_IR$res))
print(dim(phyper_IR$res))

phyper_gene <- compare_sets(df_gene_set, "gene",length(unique(df_gene_set$X)), "X")
print(head(phyper_gene$res))
print(dim(phyper_gene$res))

# Make combined cowplot with all 3 types of sets
plot_grid(phyper_gene$hist, phyper_splice$hist, phyper_IR$hist, 
        # labels = c('gene','splice', 'IR'),
        ncol = 1)

ggsave(paste0(fig_output_path,"gene_splice_IR_phyper_hist.png"), width=8, height=20, unit="cm")


comb <- combn(unique(unique(df_splice_set$set)), 2) %>%
    t() %>%
    as.data.frame()

print(ncol(comb))
print(nrow(comb))


stopCluster(cl)
