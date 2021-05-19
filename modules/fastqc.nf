process FASTQC {
  publishDir "${params.outdir}/fastqc_out"

  input:
  tuple val(name), path(reads)

  output:
  file "*_fastqc.{zip,html}"

  script:
  """
  fastqc -q $reads
  """
}
