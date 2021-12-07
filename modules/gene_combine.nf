process COMBINE_GENE {
  publishDir "${params.outdir}/combine_gene_out"
  echo true

  input:
  path metadata

  """
  # Import all kallisto run deseq2, make UMAPS before and after batch correction
  #combineGeneDE.R -m ${params.manifest} -o ${params.outdir}/combine_gene_out

  # Import combined deseq2 results and make reference matrix
  makeGeneRefMatrix.R -i ${params.outdir}/combine_gene_out  -o ${params.outdir}/combine_gene_out

  """
}
