#!/usr/bin/env Rscript

library(foreach)
library(doParallel)
cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

methods <- list("gsva","ssgsea","zscore","plage")
kcdf<- list("Gaussian","Poisson","none")
abs_ranking<- list("--abs_ranking", "")
mx_diff <- list("True", "False")

# j = kcdf, k = abs_ranking, l = mx_diff

foreach(i = methods) %do% {
    foreach(j = kcdf) %do% {
        foreach(k = abs_ranking) %do% {
            foreach(l = mx_diff) %do% {

                paste(i, j, k, l)

                outfile<- paste(i, j, k,"mxdiff", l, sep = "_")
                print(outfile)

                # Run version of gsva 
                cmd<- paste0("sudo docker run -v $PWD:/$PWD vacation/gsva:1.0.4 GSVA ",
                "--gmt /mnt/results_sigil_combine/sigil_results_SongChoi_newlabels_20220610/combine_mesa_out/splice_set_PS/gmt/main_set.gmt ",
                "--tsv_in  /mnt/results/sigil_results_SRP125125_Monaco_20220524/mesa_out/mesa_allPS_log2.tsv ",
                "--output  /mnt/monaco_benchmarking/gsva_log_outputs/",outfile,".csv", 
                " --method ", i , " --kcdf ", j , "  ", k, "  ", "--mx_diff ", l,
                " --verbose ")

                cat("\n")
                cat(cmd)
                cat("\n")
                # system(cmd)

                # make output dir for this run
                system(paste0("mkdir /mnt/monaco_benchmarking/gsva_log_eval/",outfile,"/"))

                # run script to calculate sensitivity
                cmd_eval <- paste0("sudo docker run -v $PWD:/$PWD althornt/sigl_prep:latest /mnt/sigil/bin/evalEnrichment_monaco.R ",
                    "-m /mnt/sra-manifest/SRP125125_SraRunTable_Monaco_sigil.csv ",
                    "-i /mnt/monaco_benchmarking/gsva_log_outputs/",outfile,".csv ",
                    "-o /mnt/monaco_benchmarking/gsva_log_eval/",outfile,"/")


                cat("\n")
                cat(cmd_eval)
                cat("\n")
                system(cmd_eval)



            }
        }
    }
}




parallel::stopCluster(cl)