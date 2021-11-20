#!/usr/bin/env nextflow

// enable modules
nextflow.enable.dsl = 2

// import modules
// include { KALLISTO_SE } from './modules/kallisto_se'
// include { KALLISTO_PE } from './modules/kallisto_pe'
include { KALLISTO_SE; KALLISTO_PE } from './modules/kallisto'
include { STAR_ALIGN } from './modules/star_align'
include { MESA; MESA_ONLY } from './modules/mesa'
// include { MESA_ONLY } from './modules/mesa_only'
// include { FASTQC } from './modules/fastqc'
// include { MULTIQC } from './modules/multiqc'
include { FASTQC; MULTIQC } from './modules/qc'
include { POST_KALLISTO; POST_KALLISTO_ONLY } from './modules/post_kallisto'
include { POST_MESA; POST_MESA_ONLY } from './modules/post_mesa'

// kallisto
// star
// mesa
//  QC
// post


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf

        --help                         This usage statement.
        --outdir                       Output directory to place outputs
        --reads                        Input directory of fastq.gz files
        --single_end                   Indicate reads are single end for kallisto
        --paired_end                   Indicate reads are paired end for kallis
        --kallisto_idx                 Kallisto index
        --star_index                   STAR index
        --gtf                          gtf file needed for MESA
        --metadata                     csv file needed for MESA (usually SRA run table)
        --skip_QC                      Dont run fastQC or multiQC
        --star_bed_dir                 Provide the directory with .bed files to run MESA
        --cluster                      Just run kallisto cluster step, requires kallisto directory and MESA PS

        """
}

// params.outdir = null
// params.reads = null
// params.paired_end = null
// params.single_end = null
// params.kallisto_idx = null
// params.star_index = null
// params.metadata = null
// params.skip_QC = null
// params.star_bed_dir = null

// SRP125125 (Monaco et al)
// params.outdir = "/mnt/sigil/sigil_results_20210511"
// params.reads = "/mnt/sra-fastq/*{1,2}.fastq.gz"
// params.sraruntable  = '/mnt/sra-manifest/SRP125125_SraRunTable_noblanks.csv'

// SRP253519 (Song et al)
// params.outdir = "/mnt/sigil/sigil_results_SRP253519_20210519"
// params.reads = "/mnt/sra-fastq-SRP253519/*.fastq"
// params.metadata  = "/mnt/sra-manifest/SRP253519_SraRunTable.csv"
// params.star_bed_dir = "/mnt/sigil/sigil_results_SRP253519_20210519/star_out"


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


// params.star_bed_dir = "/mnt/sigil/sigil_results_test_SE/star_out"

// params for building indices
// params.transcriptome = "/mnt/files/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.genome = "/mnt/files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
// params.annotation    = "/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf"
// params.overhang      = '99'

// using built indices
params.star_index = "/mnt/files/STARgenome"
// genome_index = file(params.star_index)
params.kallisto_idx  = '/mnt/files/homo_sapiens/transcriptome.idx'
// transcriptome_index  = file(params.kallisto_idx)

// needed for post mesa - mesa runCompareSampleSets_1_vs_all
params.gtf  = '/mnt/files/homo_sapiens/Homo_sapiens.GRCh38.96.gtf'

println "\n"
println "Input fastq directory: $params.reads "
println "Output directory: $params.outdir"
println "params.paired_end : $params.paired_end"
println "params.single_end : $params.single_end \n"



workflow {
  main:

  //kallisto PAIRED END
  if (params.paired_end && params.kallisto_idx){
    dir =  params.reads +"*{1,2}.fastq.gz"
    read_ch = channel.fromFilePairs(dir, checkIfExists: true ).view()
    KALLISTO_PE(params.kallisto_idx, read_ch)
    POST_KALLISTO(params.metadata, KALLISTO_PE.out.collect())
  }

  //kallisto SINGLE END
  if (params.single_end && params.kallisto_idx) {
    dir =  params.reads +"*.fastq.gz"
    read_ch = channel.fromPath( dir, checkIfExists: true ).map { it -> [it.name.replace(".fastq.gz", ""), file(it)]}.view()
    KALLISTO_SE(params.kallisto_idx, read_ch)
    POST_KALLISTO(params.metadata, KALLISTO_SE.out.collect())
  }

  // STAR and MESA
  if (params.star_index && !params.star_bed_dir && !params.cluster ){
    STAR_ALIGN(params.star_index, read_ch)
    MESA(params.metadata, STAR_ALIGN.out[1].collect(), params.gtf, params.genome)
    POST_MESA(params.metadata, MESA.out.collect(), params.gtf)
  }

  // MESA
  // if (params.metadata && !params.star_bed_dir && !params.cluster){
  //   MESA(params.metadata, STAR_ALIGN.out[1].collect())
  //   POST_MESA(params.metadata, MESA.out.collect())
  // }

  //Just MESA from STAR beds
  if (params.metadata && params.star_bed_dir && !params.cluster){
    MESA_ONLY(params.metadata, params.star_bed_dir ,params.gtf, params.genome)
  }

  //fastQC
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