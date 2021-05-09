#!/usr/bin/env nextflow

// enable modules
nextflow.enable.dsl = 2

//default parameters

params.reads = "/mnt/sra-fastq-test/*_{1,2}.fastq.gz"
params.transcriptome = "/mnt/files/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.outdir = "results"

// import modules
include { RNASEQ } from './modules/rnaseq'
include { STAR } from './modules/star'

workflow {
  read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
  RNASEQ( params.transcriptome, read_pairs_ch )
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
