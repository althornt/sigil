params.outdir = 'results'

include { INDEX } from './index'
include { QUANT } from './quant'

workflow RNASEQ {
  take:
    transcriptome
    read_pairs_ch

  main:
    INDEX(transcriptome)
    QUANT(INDEX.out, read_pairs_ch)

  emit:
     QUANT.out  | collect
}
