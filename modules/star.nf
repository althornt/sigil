params.outdir = 'results_star'

include { STAR_INDEX } from './star_index'
include { STAR_ALIGN } from './star_align'

workflow STAR {
  take:
    genome
    annotation
    overhang
    read_pairs_ch

  main:
    STAR_INDEX(genome,annotation,overhang)
    STAR_ALIGN(STAR_INDEX.out, read_pairs_ch)

  emit:
     STAR_ALIGN.out  | collect
}
