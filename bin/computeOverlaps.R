#!/usr/bin/env Rscript


suppressMessages({
library(optparse)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(cowplot)
})

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

# df_metadata <- read.csv(file = paste0(sigil_out_path, "combine_mesa_out/merged_metadata.csv"), stringsAsFactors=FALSE)
# ls_group_cell_types <- unlist(unique(df_metadata[["group_label"]]))
# ls_main_cell_types <- unlist(unique(df_metadata[["main_label"]]))

compare_sets <-  function(df_sets, tag,total,col){
    ls_phyper <- list()

    df_comb <- combn(unique(df_sets$set), 2) %>%
        t() %>%
        as.data.frame()

    df_out <- foreach(i = 1:nrow(df_comb),.combine='rbind', .packages=c('magrittr','dplyr','ggplot2'))  %dopar% {

        setA <- df_comb[i,"V1"]
        setB <- df_comb[i,"V2"]

        # print(paste(setA,setB,sep=" : "))
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
        g <- c(setA,setB,pval,overlap)
        }

    colnames(df_out) <- c("setA","setB","pval","overlap")
    df_out <- df_out %>% as.data.frame()
    df_out$pval = as.numeric(as.character(df_out$pval))
    df_out$overlap = as.numeric(as.character(df_out$overlap))

    df_out<- df_out %>%
        mutate(sig = ifelse(pval < .05, "significant overlap", "no significant overlap"))
    
    # Make histogram of pvalues
    hist_pval <- ggplot(df_out, aes(x=pval)) + 
        geom_histogram() +
        ggtitle(tag) +
        ylim(0, 4000) +
        theme_classic() 
    # ggsave(paste0(fig_output_path,"phyper_hist_",tag,".png"), width=10, height=5, unit="cm")

    # Make histogram of overlap
    hist_overlap <- ggplot(df_out, aes(x=overlap)) + 
        geom_histogram() +
        ggtitle(tag) +
        ylim(0, 4000) +
        theme_classic() 

    return(list("res" = df_out, "hist_pval"=hist_pval,"hist_overlap"=hist_overlap ))
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

# Make combined pval cowplot with all 3 types of sets
plot_grid(phyper_gene$hist_pval, phyper_splice$hist_pval, phyper_IR$hist_pval, 
        # labels = c('gene','splice', 'IR'),
        ncol = 1)

ggsave(paste0(fig_output_path,"gene_splice_IR_phyper_hist.png"), width=8, height=20, unit="cm")

# Make combined overlap cowplot with all 3 types of sets
plot_grid(phyper_gene$hist_overlap, phyper_splice$hist_overlap, phyper_IR$hist_overlap, 
        # labels = c('gene','splice', 'IR'),
        ncol = 1)

ggsave(paste0(fig_output_path,"gene_splice_IR_phyper_overlap.png"), width=8, height=20, unit="cm")


# Combine 3 dfs
phyper_gene$res$type <- "gene"
phyper_splice$res$type <- "splice"
phyper_IR$res$type <- "IR"
df_combined_phyper <- do.call("rbind", list(phyper_gene$res, phyper_splice$res, phyper_IR$res)) %>%
    as.data.frame()

print(head(df_combined_phyper))
ggplot(df_combined_phyper, aes(sig, ..count.., fill = type)) + 
        geom_bar( position = "dodge") +
        theme_classic() +
        theme(axis.title.x = element_blank()) +
        scale_color_manual("legend",values=c("gene"="orange","splice"="skyblue","IR"="darkgreen"))
  

ggsave(paste0(fig_output_path,"gene_splice_IR_sig.png"), width=15, height=6, unit="cm")


stopCluster(cl)


