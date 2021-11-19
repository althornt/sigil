process COMBINE_GENE {
  publishDir "${params.outdir}/combine_gene_out"
  echo true

  input:
  path metadata

  """
  combineGene.R -m ${params.manifest} -o ${params.outdir}/combine_gene_out
  """
}
