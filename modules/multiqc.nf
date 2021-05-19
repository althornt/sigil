process MULTIQC {

    publishDir "${params.outdir}/fastqc_out"

    input:
    path('*') 

    output:
    file "multiqc_report.html"

    script:
    """
    multiqc .
    """
}
