process MESA {
  echo true
  publishDir "${params.outdir}/mesa_out", mode: 'copy'

  input:
  path srafile
  file bam // not directly used, just needed so MESA runs after STAR

  output:
  //path "${params.outdir}/mesa_out"
  //file mesa_manifest.txt
  //path "${params.outdir}/mesa_out"
  //publishDir "${params.outdir}/mesa_out"
  //path "${params.outdir}/mesa_out/", emit: dir_mesa_out


  // for some reason mesa_allPS.tsv path below doesnt work but including the test file does
  //path 'mesa_allPS.tsv'
  path 'test.tsv'


  script:

  //use the srafile/metadata to create a manifest needed for mesa
  //run MESA
  //move mesa outputs to the output dir

    """
    mkdir -p ${params.outdir}/mesa_out

    sed 's/\r//' $srafile  |
    awk -F, '{print \$1",${params.outdir}/star_out/STAR_"\$1"/"\$1"Aligned.sortedByCoord.out.bam,"\$NF","\$NF}' OFS=  |     #rearrange columns, add file path
    sed '1d' |                                                                                             #remove first line
    sed -E 's/("([^"]*)")?,/\2\t/g'|                                                                       #csv to tab
    tr -cd '\11\12\15\40-\176' >  mesa_manifest.txt

    mesa quant -m mesa_manifest.txt -o mesa --drim --maxLength  50000 --minLength 50 --minOverhang 5 --minUnique 5 --lowCoverageNan --minEntropy 1

    cp mesa* ${params.outdir}/mesa_out

    echo "test" > test.tsv

    """
}
