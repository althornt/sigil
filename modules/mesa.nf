// Each process has different combinatons of MESA steps to make it easy to rerun
// individual subsets. The MESA bam_to_junc_bed command takes the longest and is
// not includes in any of the "only" processes. All commands and parameters
// should be identical.

// MESA: Runs MESA directly after STAR alignment
// MESA_ONLY: Runs MESA in a new run with existing bams from STAR
// MESA_QUANT_ONLY: Runs only MESA quant step in a new run with existing bams from STAR
// MESA_IR_ONLY: Runs only MESA quant intron retention analysis step in a new run with existing bams from STAR


process MESA {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  cpus 15

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
    # awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,"\$NF","\$NF}' OFS=  |     #rearrange columns, add file path
    awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,NA,NA"}' OFS=  |     #rearrange columns, add metadata columns as NA, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest_bams.txt

    cp mesa_manifest_bams.txt ${params.outdir}/mesa_out

    # Run MESA bam_to_junc_bed
    mesa bam_to_junc_bed -m mesa_manifest_bams.txt --output_prefix mesa \
        --annotation $gtf --genome $genome --number_threads $task.cpus \
        --strand inferCombine

    # Run MESA quant
    mesa quant -m mesa_manifest.txt -o mesa --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    # Run MESA intron coverage
    mkdir -p mesa_intron_coverage
    mesa intron_coverage -b mesa_manifest_bams.txt -m mesa_allPS.tsv \
        -j mesa_junctions.bed -n $task.cpus -o mesa_intron_coverage

    mkdir -p ${params.outdir}/mesa_out/mesa_intron_coverage
    cp *_intron_coverage.txt ${params.outdir}/mesa_out/mesa_intron_coverage

    mesa ir_table -i mesa_inclusionCounts.tsv -c mesa_allClusters.tsv \
        -d ${params.outdir}/mesa_out/mesa_intron_coverage -o mesa_ir_table -r

    cp -r mesa* ${params.outdir}/mesa_out

    """
}

process MESA_ONLY {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  cpus 15

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
        --annotation $gtf --genome $genome --number_threads $task.cpus \
        --strand inferCombine


    # Run MESA quant
    mesa quant -m mesa_manifest.txt -o mesa --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    # Run MESA intron coverage
    mkdir -p mesa_intron_coverage
    mesa intron_coverage -b mesa_manifest_bams.txt -m mesa_allPS.tsv \
        -j mesa_junctions.bed -n $task.cpus -o mesa_intron_coverage

    mkdir -p ${params.outdir}/mesa_out/mesa_intron_coverage
    cp *_intron_coverage.txt ${params.outdir}/mesa_out/mesa_intron_coverage

    mesa ir_table -i mesa_inclusionCounts.tsv -c mesa_allClusters.tsv \
        -d ${params.outdir}/mesa_out/mesa_intron_coverage -o mesa_ir_table -r

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
    mesa quant -m $bed_manifest -o mesa --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    rm mesa_manifest.txt #need to rm to avoid cp error
    cp -r mesa* ${params.outdir}/mesa_out
    """
}

process MESA_IR_ONLY {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  cpus 15

  input:
  path bed_manifest
  path junctions_bed
  path bam_manifest


  output:

  path 'error.txt'

  script:
  //run MESA quant, mesa_intron_coverage, ir_table
  //move mesa outputs to the output dir

    """
    mkdir -p ${params.outdir}/mesa_out

    mesa quant -m $bed_manifest -o mesa --maxLength  50000 \
        --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan \
        --minEntropy 1 2> error.txt

    # Run MESA intron coverage
    mkdir -p mesa_intron_coverage

    mesa intron_coverage -b $bam_manifest -m mesa_allPS.tsv \
        -j $junctions_bed -n $task.cpus -o mesa_intron_coverage

    mkdir -p ${params.outdir}/mesa_out/mesa_intron_coverage
    cp *_intron_coverage.txt ${params.outdir}/mesa_out/mesa_intron_coverage


    mesa ir_table -i mesa_inclusionCounts.tsv -c mesa_allClusters.tsv \
        -d ${params.outdir}/mesa_out/mesa_intron_coverage -o mesa_ir_table -r

    #need to rm to avoid cp error
    rm mesa_junctions.bed mesa_manifest.txt mesa_manifest_bams.txt
    cp -r mesa* ${params.outdir}/mesa_out

    """
}
