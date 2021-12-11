process COMBINE_GENE {
  publishDir "${params.outdir}/combine_gene_out"
  echo true

  input:
  path metadata

  """
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
  # combine MESA inclusion count tables


  echo "hihihihi"

  # run 1 vs all with drim seq batch correction


  # batch correct counts
  ## combine MESA counts ?
  ## log transform
  ## limma
  ## un log
  ## run MESA quant


  # Import combined deseq2 results and make reference matrix
  #makeSpliceRefMatrix.R -i ${params.outdir}/combine_mesa_out  -o ${params.outdir}/combine_mesa_out

  """
}
