# sigil
The sigil nextflow workflow uses the Docker container:  [althornt/sigl_prep](https://hub.docker.com/r/althornt/sigl_prep)

Sigil expects the metadata file to have a columns:
- "sigil_general" (more general cell type groups)
- "cell_type" (specific, treatment info)
- "LM22" (cibersort LM22 mapping)

sigil_process.nf is run per dataset

sigil_combine.nf merges all data sets

# SRA download

`bash bin/fasterqDumpFromMeta.sh /mnt/sra-manifest/SraRunTable.csv /mnt/fastq/sra-fastq 8`

# sigil_process.nf
`nextflow run sigil_process.nf`

## sigil process test data

_**Paired end**_


**Testing paired end data (full)**
```
sudo nextflow run sigil_process.nf  \
    --outdir /mnt/results/sigil_results_test_PE_new  \
    --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv  \
    --reads /mnt/fastq/sra-fastq-test-small/  \
    --paired_end
```

**Testing paired end data (MESA only)**

```
```


**Testing paired end data (cluster only)**
```
sudo nextflow run sigil_process.nf \
    --outdir /mnt/results/sigil_results_test_PE  \
    --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_mini3.csv \
    --cluster
```


_**Single end**_

**Testing single end data small files (full)**
```
sudo nextflow run sigil_process.nf \
    --outdir /mnt/results/sigil_results_test_SE_new  \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_mini3.csv \
    --reads /mnt/fastq/sra-fastq-SRP253519-mini/   \
    --single_end
```

**Testing single end data (MESA only)**

```
sudo nextflow run sigil_process.nf \
    --outdir /mnt/results/sigil_results_test_SE_newmesa  \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_mini3.csv \
    --star_bed_dir /mnt/results/sigil_results_test_SE_new/star_out \
    --skip_QC
```

**Testing single end data (cluster only)**
```
sudo nextflow run sigil_process.nf     \
  --outdir /mnt/results/sigil_results_test_SE      \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_mini3.csv     \
  --cluster`
```

## sigil process real data

### [Monaco et al](https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf)
**Running Monaco et al (full)**
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP125125_20211028 \
  --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP125125/ \
  --paired_end
```

```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP125125_Monaco_20211109 \
  --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP125125/ \
  --paired_end
```

**Running Monaco et al (MESA only)**
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP125125_Monaco_20211109 \
  --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_sigil.csv \
  --star_bed_dir /mnt/results/sigil_results_SRP125125_Monaco_20211109/star_out
  --skip_QC
```
```# new
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP125125_Monaco_20211109 \
  --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_sigil.csv \
  --star_bed_dir /mnt/results/sigil_results_SRP125125_Monaco_20211109/star_out  \
  --skip_QC
```

**Running Monaco et al (cluster only)**
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP125125_20211028 \
  --metadata  /mnt/sra-manifest/SRP125125_SraRunTable_sigil.csv \
  --cluster
```

### [Song et al](https://www.sciencedirect.com/science/article/pii/S2211124719300592?via%3Dihub)

**Running Song et al (full)**

```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP253519_20210527 \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP253519/ \
  --single_end
```

```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP253519_Song_20211114 \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP253519/ \
  --single_end
```

**Running Song et al (MESA only)**

```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP253519_20210527 \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
  --star_bed_dir /mnt/results/sigil_results_SRP253519_20210527/star_out \
  --skip_QC
```
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP253519_Song_20211114 \
  --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
  --star_bed_dir /mnt/results/sigil_results_SRP253519_Song_20211114/star_out \
  --skip_QC

```

**Running Song et al (cluster only)**
```
sudo nextflow run sigil_process.nf  \
    --outdir /mnt/results/sigil_results_SRP253519_20210527 \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
    --cluster
```

```
sudo nextflow run sigil_process.nf  \
    --outdir /mnt/results/sigil_results_SRP253519_20211028/ \
    --metadata  /mnt/sra-manifest/SRP253519_SraRunTable_sigil.csv \
    --cluster
```


### [Choi et al](https://login.ezproxy.u-pec.fr/login?qurl=https://pubmed.ncbi.nlm.nih.gov%2f30395284%2f)
**Running Choi et al (full)**
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP150419_Choi_20211124 \
  --metadata  /mnt/sra-manifest/SRP150419_SraRunTable_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP150419/ \
  --single_end
```


**Running Choi et al (MESA only)**
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP150419_Choi_20211124 \
  --metadata  /mnt/sra-manifest/SRP150419_SraRunTable_sigil.csv \
  --mesa_only
  --skip_QC
```
**Running Choi et al (cluster only)**
```
sudo nextflow run sigil_process.nf  \
    --outdir /mnt/results/sigil_results_SRP150419_Choi_20211124 \
    --metadata  /mnt/sra-manifest/SRP150419_SraRunTable_sigil.csv \
    --cluster
```


### [Calderon et al](link)

SRP156452/GSE118165 - Calderon

```
bash bin/fasterqDumpFromMeta.sh \
  /mnt/sra-manifest/SRP156452_SraRunTable.csv \
  /mnt/fastq/sra-fastq-SRP156452 \
  10
```

Run full
```
sudo nextflow run sigil_process.nf \
  --outdir /mnt/results/sigil_results_SRP156452_Calderon_20211211 \
  --metadata  /mnt/sra-manifest/SraRunTable_SRP156452_Calderon_sigil.csv \
  --reads /mnt/fastq/sra-fastq-SRP156452/ \
  --paired_end
```

_____________________________________________________________________________________________

## sigil_combine.nf
```
sudo nextflow run sigil_combine.nf \
  --manifest /mnt/files/sigil_res_manifest.txt \
  --outdir /mnt/results_sigil_combine/sigil_results_20211202

sudo nextflow run sigil_combine.nf \
  --manifest /mnt/files/sigil_res_manifest.txt \
  --outdir /mnt/results_sigil_combine/sigil_results_20211208

sudo nextflow run sigil_combine.nf \
  --manifest /mnt/files/sigil_res_manifest.txt \
  --outdir /mnt/results_sigil_combine/sigil_results_dropped_samples_20211210

sudo nextflow run sigil_combine.nf \
    --manifest /mnt/files/sigil_res_manifest.txt \
    --outdir /mnt/results_sigil_combine/sigil_results_dropped_samples_v2_20211210

# used wrong and non updated manifet??!??!
sudo nextflow run sigil_combine.nf \
    --manifest /mnt/files/sigil_res_manifest.txt \
    --outdir /mnt/results_sigil_combine/sigil_results_dropped_samples_splicing_updatedSongMonaco_20220218

sudo nextflow run sigil_combine.nf \
    --manifest /mnt/files/sigil_res_manifest_20220226.txt \
    --outdir /mnt/results_sigil_combine/sigil_results_dropped_samples_splicing_updatedSongMonaco_inferCombine_20220226

sudo nextflow run sigil_build.nf \
    --manifest /mnt/files/sigil_res_manifest_20220421.txt \
    --outdir /mnt/results_sigil_combine/sigil_results_with_Calderon_20220421


sudo nextflow run sigil_build.nf \
    --manifest /mnt/files/sigil_res_manifest_20220427_droppedsongmono.txt \
    --outdir /mnt/results_sigil_combine/sigil_results_dropped_song_mono_20220426


sudo nextflow run sigil_build.nf \
    --manifest /mnt/files/sigil_res_manifest_noSongmono_withCalderon_20220430.txt \
    --outdir /mnt/results_sigil_combine/sigil_results_noSongmono_withCalderon_20220430


sudo nextflow run sigil_build.nf  \
  --manifest /mnt/files/sigil_res_manifest_noSongmono_Calderondroppedtreatment_20220503.txt  \
  --outdir /mnt/results_sigil_combine/sigil_results_noSongmono_Calderondroppedtreatment_20220503

sudo nextflow run sigil_build.nf  \
  --manifest /mnt/files/sigil_res_manifest_newformat_20220504.txt  \
  --outdir /mnt/results_sigil_combine/sigil_results_newformat_20220504


sudo nextflow run sigil_build.nf  \
  --manifest /mnt/files/sigil_res_manifest_noSongmono_Calderontreatment_20220505.txt  \
  --outdir /mnt/results_sigil_combine/sigil_results_Calderontreatment_20220505

sudo nextflow run sigil_build.nf  \
  --manifest /mnt/sigil_build_manifests/sigil_res_grant_SongChoi_20220507.txt  \
  --outdir /mnt/results_sigil_combine/sigil_results_sigil_res_grant_SongChoi_20220507

sudo nextflow run sigil_build.nf  \
  --manifest /mnt/sigil_build_manifests/song_only.txt  \
  --outdir /mnt/results_sigil_combine/test_new_bc_song_20220508

sudo nextflow run sigil_build.nf  \
  --manifest /mnt/sigil_build_manifests/sigil_res_SongChoiMonaco_20220515.txt  \
  --outdir /mnt/results_sigil_combine/sigil_results_SongChoiMonaco_202205015_newmethod

```

# Bam to bigwigs 
`sudo docker run -v /mnt:/mnt_ -ti althornt/sigl_prep:latest` 

Song
```
./bam2bigwigs.R --manifest /mnt_/sra-manifest/SRP253519_SraRunTable_Song_grant.csv --bam_dir /mnt_/results/sigil_results_SRP253519_Song_20220507/star_out/ --out_dir /mnt_/sigil_tracks/Song/ --header /mnt_/sigil/files/header_edited.txt
```

Choi
```
./bam2bigwigs.R --manifest /mnt_/sra-manifest/SRP150419_SraRunTable_Choi_sigil.csv --bam_dir /mnt_/results/sigil_results_SRP150419_Choi_20211124/star_out --out_dir /mnt_/sigil_tracks/Choi/ --header /mnt_/sigil/files/header_edited.txt
```

Monaco
```
./bam2bigwigs.R --manifest /mnt_/sra-manifest/SRP125125_SraRunTable_Monaco_sigil.csv --bam_dir /mnt_/results/sigil_results_SRP125125_Monaco_20211109/star_out --out_dir /mnt_/sigil_tracks/Monaco/ --header /mnt_/sigil/files/header_edited.txt
```

# Make browser trackDb file 
`sudo docker run -v /mnt:/mnt_ -ti althornt/sigl_prep:latest` 

Choi
```bash
./maketrackDb.R --manifest /mnt_/sra-manifest/SRP150419_SraRunTable_Choi_sigil.csv \
     --out_dir /mnt_/sigil_tracks/Choi/ \
     --public_path "http://public.gi.ucsc.edu/~althornt/Choi_PMID30395284_tracks/"
```

Monaco 
```bash
./maketrackDb.R --manifest /mnt_/sra-manifest/SRP125125_SraRunTable_Monaco_sigil.csv \
     --out_dir /mnt_/sigil_tracks/Monaco/ \
     --public_path "http://public.gi.ucsc.edu/~althornt/Monaco_PMID30726743_tracks/"
```