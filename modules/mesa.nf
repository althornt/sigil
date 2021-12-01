process MESA {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  cpus 5

  input:
  path srafile
  file bam // not directly used, just needed so MESA runs after STAR
  path gtf
  path genome

  output:
  // for some reason using the mesa_allPS.tsv path below doesnt work
  // but including the error file does
  //path 'mesa_allPS.tsv'
  path 'error.txt'

  script:
  //use the srafile/metadata to create a manifest needed for mesa
  //run MESA
  //move mesa outputs to the output dir

    """
    mkdir -p ${params.outdir}/mesa_out

    # Make MESA manifest of bams
    sed 's/\r//' $srafile  |
    awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,"\$NF","\$NF}' OFS=  |     #rearrange columns, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest_bams.txt

    cp mesa_manifest_bams.txt ${params.outdir}/mesa_out

    # Run MESA bam_to_junc_bed
    mesa bam_to_junc_bed -m mesa_manifest_bams.txt --output_prefix mesa \
        --annotation $gtf --genome $genome --number_threads $task.cpus

    # Run MESA quant
    mesa quant -m mesa_manifest.txt -o mesa --drim --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    cp -r mesa* ${params.outdir}/mesa_out

    """
}

process MESA_ONLY {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  cpus 5

  input:
  path srafile
  path gtf
  path genome

  output:
  // for some reason using the mesa_allPS.tsv path below doesnt work
  // but including the error file does
  //path 'mesa_allPS.tsv'
  path 'error.txt'

  script:
  //use the srafile/metadata to create a manifest needed for mesa
  //run MESA
  //move mesa outputs to the output dir

    """
    mkdir -p ${params.outdir}/mesa_out

    sed 's/\r//' $srafile | head


    # Make MESA manifest of bams
    sed 's/\r//' $srafile  |
    #awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,"\$NF","\$NF}' OFS=  |     #rearrange columns, add metadata columns,  add file path
    awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,NA,NA"}' OFS=  |     #rearrange columns, add metadata columns as NA, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest_bams.txt

    cp mesa_manifest_bams.txt ${params.outdir}/mesa_out

    # Run MESA bam_to_junc_bed
    mesa bam_to_junc_bed -m mesa_manifest_bams.txt --output_prefix mesa \
        --annotation $gtf --genome $genome --number_threads $task.cpus

    # Run MESA quant
    mesa quant -m mesa_manifest.txt -o mesa --drim --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    cp -r mesa* ${params.outdir}/mesa_out

    """
}

process MESA_QUANT_ONLY {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  input:
  path bed_manifest

  output:
  path 'error.txt'

  //Run MESA quant step
    """
    mkdir -p ${params.outdir}/mesa_out

    # Run MESA quant
    mesa quant -m $bed_manifest -o mesa --drim --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    cp -r mesa* ${params.outdir}/mesa_out
    """
}
