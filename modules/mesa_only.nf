process MESA_ONLY {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  input:
  path metadata
  path star_bed_dir
  path gtf
  path genome

  output:
  path 'error.txt'

  //create manifest and run mesa bam to junction and then mesa quant
    """
    mkdir -p ${params.outdir}/mesa_out

    # Make MESA manifest of bams
    sed 's/\r//' $metadata  |
    awk -F, '{print \$1",$star_bed_dir/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,"\$NF","\$NF}' OFS=  |     #rearrange columns, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest_bams.txt

    cp mesa_manifest_bams.txt ${params.outdir}/mesa_out

    # Run MESA bam_to_junc_bed
    mesa bam_to_junc_bed -m mesa_manifest_bams.txt --output_prefix mesa --annotation $gtf --genome $genome

    #cp mesa_manifest.txt  ${params.outdir}/mesa_out
    #cp -r mesa_junction_beds/ ${params.outdir}/mesa_out

    # Run MESA quant
    mesa quant -m mesa_manifest.txt -o mesa --drim --maxLength  50000 --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan --minEntropy 1 2> error.txt
    cp -r mesa* ${params.outdir}/mesa_out


    """
}
