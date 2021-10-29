
process KALLISTO_PE {
    tag "$pair_id"
    publishDir "${params.outdir}/kallisto_out", mode: 'copy'
    cpus 5

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

    mkdir -p ${params.outdir}/kallisto_out
    cp -r $pair_id ${params.outdir}/kallisto_out
    """
}
