process FASTQC {
  publishDir "${params.outdir}/fastqc_out", mode: 'copy'

  input:
  tuple val(name), path(reads)

  output:
  file "*_fastqc.{zip,html}"

  script:
  """
  fastqc -q $reads
  """
}

process MULTIQC {

  publishDir "${params.outdir}/fastqc_out"

  input:
  path('*')

  output:
  file "multiqc_report.html"

  script:
  """
  multiqc .
  """
}
