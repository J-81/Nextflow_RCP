/*
 * Different Gene Expression Analysis Processes
 */
process DGE_BY_DESEQ2 {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${meta.DESeq2_NormCount}",
    mode: params.publish_dir_mode,
    pattern: "norm_counts_output/*", saveAs: { "${file(it).getName()}" }
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${meta.DESeq2_DGE}",
    mode: params.publish_dir_mode,
    pattern: "dge_output/*", saveAs: { "${file(it).getName()}" }
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${meta.DESeq2_DGE}/ERCC_NormDGE",
    mode: params.publish_dir_mode,
    pattern: "dge_output_ercc/*", saveAs: { "${file(it).getName()}" }

  input:
    path("runsheet.csv")
    path(organisms_csv)
    path("Rsem_gene_counts/*")
    val(meta)

  output:
    tuple path("norm_counts_output/Normalized_Counts.csv"),
          path("norm_counts_output/SampleTable.csv"),
          path("norm_counts_output/Unnormalized_Counts.csv"), emit: norm_counts

    tuple path("dge_output/contrasts.csv"),
          path("dge_output/differential_expression.csv"),
          path("dge_output/visualization_output_table.csv"),
          path("dge_output/visualization_PCA_table.csv"), emit: dge

    tuple path("norm_counts_output/ERCC_Normalized_Counts.csv"),
          path("dge_output_ercc/ERCCnorm_contrasts.csv"),
          path("dge_output_ercc/ERCCnorm_differential_expression.csv"),
          path("dge_output_ercc/visualization_output_table_ERCCnorm.csv"),
          path("dge_output_ercc/visualization_PCA_table_ERCCnorm.csv"), optional: true, emit: dge_ercc

    path("versions.txt"), emit: version

  stub:
    """
    # create output directories
    mkdir norm_counts_output
    mkdir dge_output
    ${ meta.has_ercc ? 'mkdir dge_output_ercc' : '' } # create directory, only if meta.has_ercc is true

    # run the script with R
    deseq2_normcounts_DGE_vis_ISA.R \
      ${ meta.organism_non_sci } \
      runsheet.csv \
      norm_counts_output \
      dge_output \
      TRUE \
      ${ meta.has_ercc ? 'dge_output_ercc' : '' }
    """

  script:
    """
    # create output directories
    mkdir norm_counts_output
    mkdir dge_output
    ${ meta.has_ercc ? 'mkdir dge_output_ercc' : '' } # create directory, only if meta.has_ercc is true

    # run the script with R
    deseq2_normcounts_DGE_vis_ISA.R \
      ${ meta.organism_non_sci } \
      runsheet.csv \
      norm_counts_output \
      dge_output \
      FALSE \
      ${ meta.has_ercc ? 'dge_output_ercc' : '' }
    """
}
