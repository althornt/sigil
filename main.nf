#!/usr/bin/env nextflow

// enable modules
nextflow.enable.dsl = 2

// import modules
include { KALLISTO_SE } from './modules/kallisto_se'
include { KALLISTO_PE } from './modules/kallisto_pe'
include { STAR_ALIGN_PE } from './modules/star_align_pe'
include { MESA } from './modules/mesa'


// SRP125125
// params.outdir = "/mnt/sigil/sigil_results_20210511"
// params.reads = "/mnt/sra-fastq/*{1,2}.fastq.gz"
// params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_noblanks.csv'


// pair end test files
// params.outdir = "/mnt/sigil/sigil_results_test"
// params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv'
// params.reads = "/mnt/sra-fastq-test-small/*{1,2}_short.fastq.gz"
// params.reads = "/mnt/sra-fastq-test/*{1,2}.fastq.gz"


// single end test files
params.outdir = "/mnt/sigil/sigil_results_test_SE"
params.sraruntable  = 'sra-manifest/SRP253519_SraRunTable_mini3.csv'
params.reads = "/mnt/sra-fastq-SRP253519/*.fastq"



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

params.PAIRED_END = false
params.SINGLE_END = true


workflow {
  main:

  //need indices
  if (params.star_index)
  if (params.kallisto_idx)

  //PAIRED END
  if (params.PAIRED_END){

  //kallisto and star
  read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ).view()
  KALLISTO_PE(params.kallisto_idx, read_pairs_ch)
  STAR_ALIGN(params.star_index, read_pairs_ch)

  //mesa
  if (params.sraruntable)
  MESA(params.sraruntable,
    STAR_ALIGN.out[1].collect())
  }

  //SINGLE END
  if (params.SINGLE_END){
  //kallisto and star
  read_ch = channel.fromPath( params.reads, checkIfExists: true ).map { it -> [it.name.replace(".fastq", ""), file(it)]}.view()
  KALLISTO_SE(params.kallisto_idx, read_ch)
  }

  }


workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
