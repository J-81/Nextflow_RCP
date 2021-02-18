/*
 * Different Gene Expression Analysis Processes
 */

/*
 * Script written by NASA, parameterized by dataset
 */

// TODO: rewrite for general R script, including folder naming by params (see below)!
//     # mkdir ${ params.deseq2NormPath }
//    # mkdir ${ params.deseq2DgePath }

process DGE_BY_DESEQ2 {
  conda "${baseDir}/envs/r_deseq2.yml"
  publishDir "${params.publishDirPath}/${params.deseq2NormPath}", pattern: "norm_counts_output/*"
  publishDir "${params.publishDirPath}/${params.deseq2DgePath}", pattern: "dge_output/*"


  input:
    tuple path(ISA_zip), path(organisms_csv), path(Rsem_gene_counts)
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
    # create rsem counts directory and move them into it
    mkdir RSEM_GENE_COUNTS
    mv $Rsem_gene_counts RSEM_GENE_COUNTS

    # create metadata dir and move metadata to script expected places
    mkdir metaDir
    mv $ISA_zip metaDir

    # create output directories
    mkdir norm_counts_output
    mkdir dge_output

    # run the script with R
    GLDS-104-norm_DGE_analysis.R
    """
}
