#!/usr/bin/env nextflow

// to be ran after sigil_process has been run on each data set
// combining gene and splicing outputs of multiple data sets

// enable modules
nextflow.enable.dsl = 2

// import modules
include { BUILD_GENE; BUILD_MESA; GENE_AND_SPLICING} from './modules/build'


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
  BUILD_GENE(params.manifest)

  //MESA splicing
  // BUILD_MESA(params.manifest, params.gtf)

  //Gene and MESA splicing
  // GENE_AND_SPLICING(params.manifest)
                      

  }
