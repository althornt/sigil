
process QUANT {
    tag "$pair_id"
    publishDir "${params.outdir}/kallisto_out"

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path pair_id

    println "$params.outdir"
    println publishDir

    script:
    """
    kallisto quant -i $index --threads $task.cpus -o $pair_id ${reads}
    """
}
