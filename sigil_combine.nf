#!/usr/bin/env nextflow

// combining gene and splicing outputs of multiple data sets

// enable modules
nextflow.enable.dsl = 2

// import modules
include { COMBINE_GENE; COMBINE_MESA } from './modules/combine'


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf

        --help                         This usage statement.
        --manifest
        --gtf                          gtf file needed for MESA

        """
}

params.gtf  = '/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf'

println "\n"
println "Input manifiest:  $params.manifest "
println "Output directory: $params.outdir \n"


workflow {
  main:

  //kallisto gene expression , DE, batch correction
  // COMBINE_GENE(params.manifest)

  //MESA splicing
  COMBINE_MESA(params.manifest, params.gtf)



  }
