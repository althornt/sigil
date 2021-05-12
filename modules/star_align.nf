
process STAR_ALIGN {
    echo true

    tag "$pair_id"
    publishDir "${params.outdir}/star_out"

    cpus 5

    input:
    path star_index
    tuple val(pair_id), path(reads)

    output:
    path "STAR_${pair_id}"
    file "STAR_${pair_id}/${pair_id}Aligned.sortedByCoord.out.bam"

    script:

    """

    STAR --genomeDir $star_index \
     --readFilesIn $reads \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN $task.cpus \
     --outWigType bedGraph \
     --outFileNamePrefix $pair_id \
     --readFilesCommand zcat \
     --twopassMode Basic \

     mkdir STAR_${pair_id}
     mv ${pair_id}Aligned* STAR_${pair_id}/.
     mv ${pair_id}Signal* STAR_${pair_id}/.
     mv ${pair_id}SJ* STAR_${pair_id}/.
     mv ${pair_id}Log* STAR_${pair_id}/.

     mesa star_junc_to_bed -s STAR_${pair_id}/${pair_id}SJ.out.tab

    """
}
