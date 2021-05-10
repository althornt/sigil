
process STAR_INDEX {
    tag "$genome.simpleName"

    input:
    path genome
    path annotation
    val overhang

    output:
    path 'STARgenome'

    script:
    """
    mkdir STARgenome

    STAR --runThreadN 15 \
         --runMode genomeGenerate \
         --genomeFastaFiles $genome \
         --sjdbGTFfile $annotation \
         --sjdbOverhang $params.overhang \
         --genomeDir STARgenome \
         --outFileNamePrefix STARgenome

    """
}
