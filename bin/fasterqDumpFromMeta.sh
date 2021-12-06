#!/bin/bash

metadata_run_csv="$1"
out_dir="$2"
threads="$3"

# Make output directory
mkdir -p "$out_dir"
cd "$out_dir"

# Extract SRAs from first column and download each fastq
for sra in $(cut -d, -f1 "$metadata_run_csv" | tail -n +2); do

  echo "$sra"

  # Download fastq
  cmd="sudo docker run -t --rm \
      -v $PWD:/mnt/sigil:rw \
      -w /mnt/sigil \
      ncbi/sra-tools fasterq-dump -e "$threads" -p "$sra""

  echo "$cmd"
  $cmd

  # gzip downloaded fastq with pigz (faster than gzip)
  # it will not gzip .gz files
  pigz $sra*

done
