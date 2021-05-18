
process KALLISTO_SE {
    tag "$id"
    publishDir "${params.outdir}/kallisto_out"
    cpus 5

    input:
    path index
    tuple val(id), path(reads)


    output:
    path id

    println "$params.outdir"
    println publishDir

    script:
    """
    kallisto quant -i $index --threads $task.cpus --single -l 200 -s 20 -o $id ${reads}
    """
}
