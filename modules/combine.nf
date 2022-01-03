process COMBINE_GENE {
  publishDir "${params.outdir}/combine_gene_out"
  echo true

  input:
  path metadata

  """
  echo "Running COMBINE_GENE..."

  # Import all kallisto run deseq2, make UMAPS before and after batch correction
  combineGeneDE.R -m ${params.manifest} -o ${params.outdir}/combine_gene_out

  # Import combined deseq2 results and make reference matrix
  makeGeneRefMatrix.R -i ${params.outdir}/combine_gene_out  -o ${params.outdir}/combine_gene_out

  """
}


process COMBINE_MESA {
  publishDir "${params.outdir}/combine_mesa_out"
  echo true

  input:
  path metadata

  """
  echo "Running COMBINE_MESA..."

  mkdir ${params.outdir}/combine_mesa_out -p


  # Import and combine
  combineMESAbatchcorr.R  -m ${params.manifest} -o ${params.outdir}/combine_mesa_out

  """
}
