#!/usr/bin/env nextflow

// enable modules
nextflow.enable.dsl = 2

// import modules
include { KALLISTO_SE; KALLISTO_PE } from './modules/kallisto'
include { STAR_ALIGN } from './modules/star_align'
include { MESA; MESA_ONLY; MESA_QUANT_ONLY } from './modules/mesa'
include { FASTQC; MULTIQC } from './modules/qc'
include { POST_KALLISTO; POST_KALLISTO_ONLY } from './modules/post_kallisto'
include { POST_MESA; POST_MESA_ONLY } from './modules/post_mesa'

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run sigil_process.nf

        --outdir                       Output directory to place outputs
        --reads                        Input directory of fastq.gz files
        --single_end                   Indicate reads are single end for kallisto
        --paired_end                   Indicate reads are paired end for kallis
        --kallisto_idx                 Kallisto index
        --star_index                   STAR index
        --gtf                          gtf file needed for MESA
        --metadata                     csv file needed for MESA (usually SRA run table)
        --skip_QC                      Dont run fastQC or multiQC
        --bed_manifest                 manifest with .bed files to run MESA quant
        --cluster                      Just run kallisto cluster step, requires kallisto directory and MESA PS

        """
}

// pair end test files
// params.outdir = "/mnt/sigil/sigil_results_test_PE_4"
// params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv'
// params.reads = "/mnt/sra-fastq-test-small/*{1,2}_short.fastq.gz"
// params.reads = "/mnt/sra-fastq-test/*{1,2}.fastq.gz"

// single end test files
// params.outdir = "/mnt/sigil/sigil_results_test_SE_3"
// params.metadata  = "/mnt/sigil/sra-manifest/SRP253519_SraRunTable_mini3.csv"
// params.reads = "/mnt/sra-fastq-SRP253519-mini/*.fastq"
// params.reads = "/mnt/sra-fastq-SRP253519-mini/"

// params for building indices
// params.transcriptome = "/mnt/files/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.genome = "/mnt/files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
// params.annotation    = "/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf"
// params.overhang      = '99'

// Using built indices
params.star_index = "/mnt/files/STARgenome"
params.kallisto_idx  = '/mnt/files/homo_sapiens/transcriptome.idx'
params.gtf  = '/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf'

println "\n"
println "Input fastq directory: $params.reads "
println "Output directory: $params.outdir \n"


workflow {
  main:

  // kallisto PAIRED END
  if (params.paired_end && params.kallisto_idx){
    dir =  params.reads +"*{1,2}.fastq.gz"
    read_ch = channel.fromFilePairs(dir, checkIfExists: true ).view()
    KALLISTO_PE(params.kallisto_idx, read_ch)
    POST_KALLISTO(params.metadata, KALLISTO_PE.out.collect())
  }

  // kallisto SINGLE END
  if (params.single_end && params.kallisto_idx) {
    dir =  params.reads +"*.fastq.gz"
    read_ch = channel.fromPath( dir, checkIfExists: true ).map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}.view()
    KALLISTO_SE(params.kallisto_idx, read_ch)
    POST_KALLISTO(params.metadata, KALLISTO_SE.out.collect())
  }

  // STAR and MESA
  if (params.star_index && !params.bed_manifest && !params.cluster && !params.mesa_only ){
    STAR_ALIGN(params.star_index, read_ch)
    MESA(params.metadata, STAR_ALIGN.out[1].collect(), params.gtf, params.genome)
    POST_MESA(params.metadata, MESA.out.collect(), params.gtf)
  }

  // Just MESA to generate beds from bams and quantify
  if (params.mesa_only && params.metadata && !params.cluster){
    MESA_ONLY(params.metadata, params.gtf, params.genome)
  }

  // Just MESA QUANT from MESA generated beds
  if (params.metadata && params.bed_manifest  && !params.mesa_only && !params.cluster){
    MESA_QUANT_ONLY(params.bed_manifest)
  }

  // fastQC
  if(!params.skip_QC && params.reads && !params.cluster){
    FASTQC(read_ch)
    MULTIQC(FASTQC.out.collect())
  }

  // cluster only from MESA PS and kallisto outputs
  if(params.cluster && params.metadata){
    POST_KALLISTO_ONLY(params.metadata)
    POST_MESA_ONLY(params.metadata, params.gtf)
  }

  }
