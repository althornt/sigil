process MESA_ONLY {
  publishDir "${params.outdir}/mesa_out"

  echo true

  input:
  path metadata
  path star_bed_dir

    """
    sed 's/\r//' $metadata  |
    awk -F, '{print \$1",$star_bed_dir/STAR_"\$1"/"\$1"SJ.out.tab.bed,"\$NF","\$NF}' OFS=  |     #rearrange columns, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest.txt

    mesa quant -m mesa_manifest.txt -o mesa --drim
    mkdir -p ${params.outdir}/mesa_out
    mv mesa* ${params.outdir}/mesa_out

    """
}
