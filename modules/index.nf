
process INDEX {
    tag "$transcriptome.simpleName"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    kallisto index $transcriptome -i index
    """
}
