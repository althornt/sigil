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
  path gtf

  """
  echo "Running COMBINE_MESA..."

  mkdir ${params.outdir}/combine_mesa_out -p


  # Import,  combine, batch correct
  #combineMESAbatchcorr.R  -m ${params.manifest} -o ${params.outdir}/combine_mesa_out

  # Run MESA on batch corrected values
  # command

  # Run compare sample sets on combined all PS
  # eventually change to import batch corrected mesa_allPS
  #runMESAcompare.R -i ${params.outdir}/combine_mesa_out/LM22_mesa_allPS.tsv \
  #  -o ${params.outdir}/combine_mesa_out \
  #  -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
  #  --gtf $gtf

  # Explore LM22 results
  LM22_splicing_explore.R \
    -i ${params.outdir}/combine_mesa_out/LM22_mesa_allPS.tsv \
    -o ${params.outdir}/combine_mesa_out \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv

  # Make ref matix

  """
}
