
process STAR_ALIGN {
    echo true

    tag "$pair_id"
    publishDir "${params.outdir}/star_out"


    input:
    path star_index
    tuple val(pair_id), path(reads)

    output:
    path "STAR_${pair_id}"

    script:

    // --twopassMode Basic \

    """

    STAR --genomeDir $star_index \
     --readFilesIn $reads \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 5 \
     --outWigType bedGraph \
     --outFileNamePrefix $pair_id \
     --readFilesCommand zcat \


     mkdir STAR_${pair_id}
     mv ${pair_id}Aligned* STAR_${pair_id}/.
     mv ${pair_id}Signal* STAR_${pair_id}/.
     mv ${pair_id}SJ* STAR_${pair_id}/.
     mv ${pair_id}Log* STAR_${pair_id}/.

     ls STAR_${pair_id}
     mesa star_junc_to_bed -s STAR_${pair_id}/${pair_id}SJ.out.tab
     ls STAR_${pair_id}

     head STAR_${pair_id}/${pair_id}SJ.out.tab.bed


    """
}
