# sigil

Common commands on openstack...

`sudo nextflow run main.nf`



**Testing paired end data (full)**
`sudo nextflow run main.nf \
    --outdir /mnt/sigil/sigil_results_test_PE  \
    --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv \
    --reads /mnt/sra-fastq-test-small/ \
    --paired_end   `


**Testing paired end data (MESA only)**





Single end

**Testing single end data (full)**
`sudo nextflow run main.nf \
    --outdir /mnt/sigil/sigil_results_test_SE  \
    --metadata  /mnt/sigil/sra-manifest/SRP253519_SraRunTable_mini3.csv \
    --reads /mnt/sra-fastq-SRP253519-mini/   \
    --single_end   `


**Testing single end data (MESA only)**

`sudo nextflow run main.nf \
    --outdir /mnt/sigil/sigil_results_test_SE  \
    --metadata  /mnt/sigil/sra-manifest/SRP253519_SraRunTable_mini3.csv \
    --star_bed_dir /mnt/sigil/sigil_results_test_SE/star_out \
    --skip_QC `


Song et al

**Running Song et al (full)**


**Running Song et al (MESA only)**

`sudo nextflow run main.nf
  --outdir /mnt/sigil/sigil_results_SRP253519_20210519 \
  --metadata  /mnt/sigil/sra-manifest/SRP253519_SraRunTable.csv \
  --star_bed_dir /mnt/sigil/sigil_results_SRP253519_20210519/star_out \
  --skip_QC \`
