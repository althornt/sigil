#!/usr/bin/env nextflow

// combining gene and splicing outputs of multiple data sets

// enable modules
nextflow.enable.dsl = 2

// import modules
include { COMBINE_GENE } from './modules/gene_combine'


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf

        --help                         This usage statement.
        --manifest

        """
}


println "\n"
println "Input manifiest:  $params.manifest "
println "Output directory: $params.outdir \n"



workflow {
  main:

  //kallisto gene expression
  COMBINE_GENE(params.manifest)

  //MESA splicing



  }