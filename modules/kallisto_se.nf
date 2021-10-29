
process KALLISTO_SE {
    tag "$id"
    publishDir "${params.outdir}/kallisto_out" , mode: 'copy'
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

    mkdir -p ${params.outdir}/kallisto_out
    cp -r $id/ ${params.outdir}/kallisto_out
    """
}
