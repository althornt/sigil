process POST_MESA {
  publishDir "${params.outdir}/post_mesa_out"
  echo true

  input:
  path metadata
  file mesa
  path gtf

  """
  postMESA.R -i ${params.outdir}/mesa_out/mesa_allPS.tsv -o ${params.outdir}/post_mesa_out -m $metadata
  runMESAcompare.R -i ${params.outdir}/mesa_out/mesa_allPS.tsv -o ${params.outdir}/post_mesa_out -m $metadata --gtf $gtf
  """
}


process POST_MESA_ONLY {
  publishDir "${params.outdir}/post_mesa_out"
  echo true

  input:
  path metadata
  path gtf

  """
  postMESA.R -i ${params.outdir}/mesa_out/mesa_allPS.tsv -o ${params.outdir}/post_mesa_out -m $metadata
  runMESAcompare.R -i ${params.outdir}/mesa_out/mesa_allPS.tsv -o ${params.outdir}/post_mesa_out -m $metadata --gtf $gtf
  """
}
