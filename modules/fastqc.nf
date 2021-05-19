
process FASTQC {

  tag "$id"
  publishDir "${params.outdir}/fastqc_out"

  input:
  tuple val(name), file(reads)

  output:
  file "*_fastqc.{zip,html}" 

  script:
  """
  fastqc -q $reads
  """
}
