# sigil

`nextflow run main.nf`

The sigil nextflow workflow uses the Docker container:  [althornt/sigl_prep](https://hub.docker.com/r/althornt/sigl_prep)

sigil clustering expects sra metadata file to have a columns "general" (more general cell type groups) and "cell_type" (specific, treatment info)

# SRA download

`bash bin/fastqDumpFromMeta.sh sra_manifest.csv`





# Test Data

_**Paired end**_


**Testing paired end data (full)**
```
sudo nextflow run main.nf  \
    --outdir /mnt/results/sigil_results_test_PE  \
    --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv  \
    --reads /mnt/fastq/sra-fastq-test-small/  \
    --paired_end
```


**Testing paired end data (MESA only)**

```
```


**Testing paired end data (cluster only)**
```
sudo nextflow run main.nf \
    --outdir /mnt/results/sigil_results_test_PE  \
    --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv \
    --cluster
```


_**Single end**_

**Testing single end data small files (full)**
```
sudo nextflow run main.nf \
    --outdir /mnt/results/sigil_results_test_SE  \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_mini3.csv \
    --reads /mnt/fastq/sra-fastq-SRP253519-mini/   \
    --single_end
```

**Testing single end data (MESA only)**

```
sudo nextflow run main.nf \
    --outdir /mnt/results/sigil_results_test_SE  \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_mini3.csv \
    --star_bed_dir /mnt/results/sigil_results_test_SE/star_out \
    --skip_QC `
```

**Testing single end data (cluster only)**
```
sudo nextflow run main.nf     \
  --outdir /mnt/results/sigil_results_test_SE      \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_mini3.csv     \
  --cluster`
```

# Real Data

### [Monaco et al](https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf)
**Running Monaco et al (full)**
```

```

**Running Monaco et al (MESA only)**
```

```

**Running Monaco et al (cluster only)**
```

```


### [Song et al](https://www.sciencedirect.com/science/article/pii/S2211124719300592?via%3Dihub)

**Running Song et al (full)**

```
sudo nextflow run main.nf \
  --outdir /mnt/results/sigil_results_SRP253519_20210527 \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP253519/ \
  --single_end
```

**Running Song et al (MESA only)**

```
sudo nextflow run main.nf \
  --outdir /mnt/results/sigil_results_SRP253519_20210527 \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
  --star_bed_dir /mnt/results/sigil_results_SRP253519_20210527/star_out \
  --skip_QC
```

  **Running Song et al (cluster only)**
```
sudo nextflow run main.nf  \
    --outdir /mnt/results/sigil_results_SRP253519_20210527 \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
    --cluster
```
