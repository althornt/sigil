process COMBINE_GENE {
  publishDir "${params.outdir}/combine_gene_out"
  echo true

  input:
  path metadata

  """
  echo "Running COMBINE_GENE..."

  # Import all kallisto, run deseq2, make UMAPS before and after batch correction
  combineGeneDE.R -m ${params.manifest} -o ${params.outdir}/combine_gene_out

  # Import combined deseq2 results and make reference matrix
  makeGeneRefMatrix.R -i ${params.outdir}/combine_gene_out  -o ${params.outdir}/combine_gene_out

  """
}


process COMBINE_MESA {
  publishDir "${params.outdir}/combine_mesa_out"
  echo true

  input:
  path metadata
  path gtf

  """
  echo "Running COMBINE_MESA................................................ "

  # mkdir ${params.outdir}/combine_mesa_out -p

  # Import, combine data sets, batch correct
  echo "Running combineMESAbatchcorr.R ................................................ "
  combineMESAbatchcorr.R  -m ${params.manifest} -o ${params.outdir}/combine_mesa_out

  # Run compare sample sets on batch corrected PS values using LM22 and LM6 cell types
  echo "Running runMESAcompare.R ................................................ "
  runMESAcompare.R -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv \
    -o ${params.outdir}/combine_mesa_out \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
    --gtf $gtf

  # More specific splicing comparisons within each cell type
  echo "Running compareWithinType.R ................................................ "
  compareWithinType.R -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv \
    -o ${params.outdir}/combine_mesa_out/compare_within_type \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
    --gtf $gtf

  # Explore LM22 results and make ref matrix
  # Analyze splicing comparison outputs
  echo "Running makeSplicingRefMatrix.R ................................................ "
  makeSplicingRefMatrix.R \
    -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv \
    -o ${params.outdir}/combine_mesa_out \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv

  """
}

process GENE_AND_SPLICING {
  publishDir "${params.outdir}/combine_gene_and_mesa"
  echo true

  input:
  path metadata

  """
  echo "Running GENE_AND_SPLICING"

  # Compare gene vs splicing
  #are the dfiferntailly spliced junctions also differentially expressed ?
  # Make combined gene and splicing makeGeneRefMatrix



  """
}
