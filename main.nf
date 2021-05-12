#!/usr/bin/env nextflow

// enable modules
nextflow.enable.dsl = 2

// import modules
include { QUANT } from './modules/quant'
include { STAR_ALIGN } from './modules/star_align'
include { MESA } from './modules/mesa'

//default parameters
params.outdir = "/mnt/sigil/sigil_results_test"
// params.reads = "/mnt/sra-fastq-test-small/*{1,2}_short.fastq.gz"
// params.reads = "/mnt/sra-fastq/*{1,2}.fastq.gz"
params.reads = "/mnt/sra-fastq-test/*{1,2}.fastq.gz"
params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv'

// //building indices
// //kallisto
// params.transcriptome = "/mnt/files/Homo_sapiens.GRCh38.cdna.all.fa.gz"
// //STAR
// params.genome = "/mnt/files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
// params.annotation    = "/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf"
// params.overhang      = '99'

// not building indices
params.star_index = "/mnt/files/STARgenome"
genome_index = file(params.star_index)
params.kallisto_idx  = '/mnt/files/homo_sapiens/transcriptome.idx'
transcriptome_index  = file(params.kallisto_idx)


workflow {
  main:

  if (params.star_index)
  if (params.kallisto_idx)
  //align()
  read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ).view()
  QUANT(params.kallisto_idx, read_pairs_ch)
  STAR_ALIGN(params.star_index, read_pairs_ch)

  if (params.sraruntable)
  MESA(params.sraruntable,
    STAR_ALIGN.out[1].collect())
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
