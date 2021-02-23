/*
 * Different Gene Expression Analysis Processes
 */

/*
 * Script written by NASA, parameterized by dataset
 */

// TODO: REMOVE HARD CODED CONDA PATH

process DGE_BY_DESEQ2 {
  conda "${baseDir}/envs/RNAseq_Rtools.yml"
  publishDir "${params.publishDirPath}/${params.deseq2NormPath}", pattern: "norm_counts_output/*", saveAs: { "${file(it).getName()}" }
  publishDir "${params.publishDirPath}/${params.deseq2DgePath}", pattern: "dge_output/*", saveAs: { "${file(it).getName()}"}

  input:
    tuple path(Isa_zip), path(organisms_csv), path(Rsem_gene_counts)
  output:
    tuple path("norm_counts_output/Normalized_Counts.csv"),
          path("norm_counts_output/SampleTable.csv"),
          path("norm_counts_output/Unnormalized_Counts.csv"), emit: norm_counts
    tuple path("dge_output/contrasts.csv"),
          path("dge_output/differential_expression.csv"),
          path("dge_output/visualization_output_table.csv"),
          path("dge_output/visualization_PCA_table.csv"), emit: dge

  stub:
    """
    mkdir norm_counts_output
    touch "norm_counts_output/Normalized_Counts.csv" \
          "norm_counts_output/SampleTable.csv" \
          "norm_counts_output/Unnormalized_Counts.csv"

    mkdir dge_output
    touch "dge_output/contrasts.csv" \
          "dge_output/differential_expression.csv" \
          "dge_output/visualization_output_table.csv" \
          "dge_output/visualization_PCA_table.csv"
    """

  script:
    """
    # create output directories
    mkdir norm_counts_output
    mkdir dge_output
    mkdir dge_output_ercc

    # run the script with R
    deseq2_normcounts_wERCC_DGE_vis_ISA.R \
      ${ params.organism_nonsci } \
      $Isa_zip \
      norm_counts_output \
      dge_output \
      dge_output_ercc
    """
}
