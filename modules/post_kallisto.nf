process POST_KALLISTO {
  echo true
  publishDir "${params.outdir}/post_kallisto_out"

  input:
  path metadata
  file kallisto

  """
  postKallisto.R -i ${params.outdir}/kallisto_out -o ${params.outdir}/post_kallisto_out -m $metadata
  """
}

process POST_KALLISTO_ONLY {
  echo true
  publishDir "${params.outdir}/post_kallisto_out"

  input:
  path metadata

  """
  postKallisto.R -i ${params.outdir}/kallisto_out -o ${params.outdir}/post_kallisto_out -m $metadata
  """
}
