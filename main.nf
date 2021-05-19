#!/usr/bin/env nextflow

// enable modules
nextflow.enable.dsl = 2

// import modules
include { KALLISTO_SE } from './modules/kallisto_se'
include { KALLISTO_PE } from './modules/kallisto_pe'
include { STAR_ALIGN } from './modules/star_align'
include { MESA } from './modules/mesa'
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'


// //building indices
// //kallisto
// params.transcriptome = "/mnt/files/Homo_sapiens.GRCh38.cdna.all.fa.gz"
// //STAR
// params.genome = "/mnt/files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
// params.annotation    = "/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf"
// params.overhang      = '99'


// SRP125125
// params.outdir = "/mnt/sigil/sigil_results_20210511"
// params.reads = "/mnt/sra-fastq/*{1,2}.fastq.gz"
// params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_noblanks.csv'


// pair end test files
params.outdir = "/mnt/sigil/sigil_results_test_PE"
params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv'
params.reads = "/mnt/sra-fastq-test-small/*{1,2}_short.fastq.gz"
// params.reads = "/mnt/sra-fastq-test/*{1,2}.fastq.gz"


// single end test files
// params.outdir = "/mnt/sigil/sigil_results_test_SE"
// params.sraruntable  = "/mnt/sigil/sra-manifest/SRP253519_SraRunTable_mini3.csv"
// params.reads = "/mnt/sra-fastq-SRP253519-mini/*.fastq.gz"

params.PAIRED_END = true
params.SINGLE_END = false

// not building indices
params.star_index = "/mnt/files/STARgenome"
genome_index = file(params.star_index)
params.kallisto_idx  = '/mnt/files/homo_sapiens/transcriptome.idx'
transcriptome_index  = file(params.kallisto_idx)

workflow {
  main:

  //need indices
  if (params.star_index)
  if (params.kallisto_idx)

  //kallisto PAIRED END
  if (params.PAIRED_END){
    read_ch = channel.fromFilePairs( params.reads, checkIfExists: true ).view()
    KALLISTO_PE(params.kallisto_idx, read_ch)
  }

  //kallisto SINGLE END
  if (params.SINGLE_END){
    read_ch = channel.fromPath( params.reads, checkIfExists: true ).map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}.view()
    KALLISTO_SE(params.kallisto_idx, read_ch)
  }

  //STAR
  STAR_ALIGN(params.star_index, read_ch)

  //MESA
  if (params.sraruntable)
  MESA(params.sraruntable,
    STAR_ALIGN.out[1].collect())

  //fastQC
  FASTQC(read_ch)
  MULTIQC(FASTQC.out)


  }




workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/fastqcmultiqc_report.html\n" : "Oops .. something went wrong" )
}
