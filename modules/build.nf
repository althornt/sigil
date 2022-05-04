process BUILD_GENE {
  publishDir "${params.outdir}/combine_gene_out"
  echo true

  input:
  path metadata

  """
  echo "Running COMBINE_GENE..."
  mkdir ${params.outdir}/combine_gene_out -p

  # Import all kallisto, run deseq2, make UMAPS before and after batch correction
  combineGeneDE.R -m ${params.manifest} -o ${params.outdir}/combine_gene_out

  # Import combined deseq2 results and make reference matrix
  makeGeneRefMatrix.R -i ${params.outdir}/combine_gene_out  -o ${params.outdir}/combine_gene_out

  # Explore cell type specific genes 
  exploreGeneRefMatrix.R \
   --geneRefMatrix  ${params.outdir}/combine_gene_out/ref_matrix/lm22_lm6_withinType_combinedRefMat.tsv \
   -o ${params.outdir}/combine_gene_out/explore_ref_matrix \
   -m ${params.outdir}/combine_gene_out/metadata.csv \
   -i ${params.outdir}/combine_gene_out/combined_kallisto_log2tpm_batch_corrected.csv

  """
}


process BUILD_MESA {
  publishDir "${params.outdir}/combine_mesa_out"
  echo true

  input:
  path metadata
  path gtf

  
  // mkdir ${params.outdir}/combine_mesa_out -p

  // # Import, combine data sets, batch correct
  // combineMESAbatchcorr.R  -m ${params.manifest} -o ${params.outdir}/combine_mesa_out

  // # IR -------------------------------------------------------------------------
  // # IR Run mesa compare 
  // runMESAcompare.R -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_ir_table_intron_retention_LM22.tsv \
  //   -o ${params.outdir}/combine_mesa_out/compare_groups_IR \
  //   -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
  //   --gtf $gtf

  // # IR More specific splicing comparisons within each cell type
  // echo "Running compareWithinType.R .........................................."
  // compareWithinType.R -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_ir_table_intron_retention_LM22.tsv \
  //   -o ${params.outdir}/combine_mesa_out/compare_within_type_IR \
  //   -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
  //   --gtf $gtf
  """
  echo "Running makeIntronRetentionRefMatrix.R ......................................"
  makeIntronRetentionRefMatrix.R \
    -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_ir_table_intron_retention_LM22.tsv \
    -g ${params.outdir}/combine_mesa_out/compare_groups_IR \
    -w ${params.outdir}/combine_mesa_out/compare_within_type_IR \
    -c ${params.outdir}/combine_mesa_out/batch_corr_mesa_allClusters.tsv \
    -o ${params.outdir}/combine_mesa_out/ref_matrix_IR \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv

  
  #Explore IR ref 
  exploreSplicingRefMatrix.R \
   --spliceRefMatrix  ${params.outdir}/combine_mesa_out/ref_matrix_IR/lm22_lm6_withinType_combinedRefMat.tsv \
   -o ${params.outdir}/combine_mesa_out/explore_ref_matrix_IR \
   -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
   -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_ir_table_intron_retention_LM22.tsv

  
  # PS -------------------------------------------------------------------------

  # Run compare sample sets on batch corrected PS values using LM22 and LM6 cell types
  echo "Running runMESAcompare.R ................................................ "
  runMESAcompare.R -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv \
    -o ${params.outdir}/combine_mesa_out/compare_groups_PS \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
    --gtf $gtf

  # More specific splicing comparisons within each cell type
  echo "Running compareWithinType.R ................................................ "
  compareWithinType.R -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv \
    -o ${params.outdir}/combine_mesa_out/compare_within_type_PS \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
    --gtf $gtf

  # Explore LM22 results and make ref matrix
  # Analyze splicing comparison outputs
  echo "Running makeSplicingRefMatrix.R ................................................ "
  makeSplicingRefMatrix.R \
    -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv \
    -g ${params.outdir}/combine_mesa_out/compare_groups_PS \
    -w ${params.outdir}/combine_mesa_out/compare_within_type_PS \
    -c ${params.outdir}/combine_mesa_out/batch_corr_mesa_allClusters.tsv \
    -o ${params.outdir}/combine_mesa_out/ref_matrix_PS \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv
  
  exploreSplicingRefMatrix.R \
   --spliceRefMatrix  ${params.outdir}/combine_mesa_out/ref_matrix_PS/lm22_lm6_withinType_combinedRefMat.tsv \
   -o ${params.outdir}/combine_mesa_out/explore_ref_matrix_PS \
   -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \
   -i ${params.outdir}/combine_mesa_out/batch_corr_mesa_allPS_LM22.tsv 
  """
  }

process GENE_AND_SPLICING {
  publishDir "${params.outdir}/combine_gene_and_mesa"
  echo true

  input:
  path metadata

  """
  echo "Running GENE_AND_SPLICING"

  mkdir ${params.outdir}/gene_and_splicing_out -p

  geneAndSplicing.R \
    --spliceDir  ${params.outdir}/combine_mesa_out/ \
    --geneDir  ${params.outdir}/combine_gene_out/ \
    -o ${params.outdir}/gene_and_splicing_out/ \
    -m ${params.outdir}/combine_mesa_out/lm22_metadata.csv \


"""
}
