
process STAR_ALIGN {
    echo true

    tag "$file_id"
    publishDir "${params.outdir}/star_out" , mode: 'copy'

    cpus 5

    input:
    path star_index
    tuple val(file_id), path(reads)

    output:
    path "STAR_${file_id}"
    file "STAR_${file_id}/${file_id}Aligned.sortedByCoord.out.bam"

    script:

    """

    STAR --genomeDir $star_index \
     --readFilesIn $reads \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN $task.cpus \
     --outWigType bedGraph \
     --outFileNamePrefix $file_id \
     --twopassMode Basic \
     --readFilesCommand zcat \

     mkdir STAR_${file_id}
     mv ${file_id}Aligned* STAR_${file_id}/.
     mv ${file_id}Signal* STAR_${file_id}/.
     mv ${file_id}SJ* STAR_${file_id}/.
     mv ${file_id}Log* STAR_${file_id}/.



#     mesa star_junc_to_bed -s STAR_${file_id}/${file_id}SJ.out.tab

     # check if tab file is empty
#     if ! [ -s "STAR_${file_id}/${file_id}SJ.out.tab" ];then
#      echo "Exiting due to empty STAR .tab file"
#      exit
#     fi

     # check if bed file is empty
#     if ! [ -s "STAR_${file_id}/${file_id}SJ.out.tab.bed" ];then
#      echo "Exiting due to empty bed file from MESA star_junc_to_bed"
#      # exit
#     fi

    """
}
