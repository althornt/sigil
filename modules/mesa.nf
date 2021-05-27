process MESA {
  publishDir "${params.outdir}/mesa_out"

  echo true

  input:
  path srafile
  file bam

  output:
  path "${params.outdir}/mesa_out"

    """
    sed 's/\r//' $srafile  |
    awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"SJ.out.tab.bed,"\$NF","\$NF}' OFS=  |     #rearrange columns, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest.txt                                                        #remove weird characters, write to file

    mesa quant -m mesa_manifest.txt -o mesa --drim
    mkdir ${params.outdir}/mesa_out
    mv mesa* ${params.outdir}/mesa_out

    """
}
