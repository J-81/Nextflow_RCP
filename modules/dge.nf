/*
 * Different Gene Expression Analysis Processes
 */
process DGE_BY_DESEQ2 {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/04-DESeq2_NormCounts",
    mode: params.publish_dir_mode,
    pattern: "norm_counts_output/*", saveAs: { "${file(it).getName()}" }
  publishDir "${ params.outputDir }/${ params.gldsAccession }/05-DESeq2_DGE",
    mode: params.publish_dir_mode,
    pattern: "dge_output/*", saveAs: { "${file(it).getName()}" }
  publishDir "${ params.outputDir }/${ params.gldsAccession }/05-DESeq2_DGE",
    mode: params.publish_dir_mode,
    pattern: "dge_output_ercc/*", saveAs: { "ERCC_NormDGE/${file(it).getName()}" }

  input:
    path("runsheet.csv")
    path(organisms_csv)
    path("Rsem_gene_counts/*")
    val(meta)
    path(annotation_file)
    path("dge_annotation_R_scripts")

  output:
    tuple path("norm_counts_output/Normalized_Counts.csv"),
          path("norm_counts_output/SampleTable.csv"),
          path("norm_counts_output/RSEM_Unnormalized_Counts.csv"), emit: norm_counts

    tuple path("dge_output/contrasts.csv"),
          path("dge_output/differential_expression.csv"),
          path("dge_output/visualization_output_table.csv"),
          path("dge_output/visualization_PCA_table.csv"), emit: dge

    tuple path("norm_counts_output/ERCC_Normalized_Counts.csv"),
          path("norm_counts_output/ERCC_SampleTable.csv"),
          path("dge_output_ercc/ERCCnorm_contrasts.csv"),
          path("dge_output_ercc/ERCCnorm_differential_expression.csv"),
          path("dge_output_ercc/visualization_output_table_ERCCnorm.csv"),
          path("dge_output_ercc/visualization_PCA_table_ERCCnorm.csv"), optional: true, emit: dge_ercc

    path("versions.txt"), emit: version

  stub:
    """
    ./dge_annotation_R_scripts/dge_annotation_workflow.R \\
        --runsheet_path runsheet.csv \\
        --DEBUG_MODE_ADD_DUMMY_COUNTS \\
        --input_gene_results_dir "Rsem_gene_counts" \\
        --primary_keytype ${ meta.primary_keytype } \\
        --normalization 'default' \\
        --normalized_counts_output_prefix "norm_counts_output/" \\
        --dge_output_prefix "dge_output/" \\
        --annotation_file_path ${annotation_file} \\
        --extended_table_output_prefix "dge_output/"\\
        --extended_table_output_suffix ".csv" \\
        --verbose

    if ${ meta.has_ercc ? 'true' : 'false'}
    then
        ./dge_annotation_R_scripts/dge_annotation_workflow.R \\
            --runsheet_path runsheet.csv \\
            --DEBUG_MODE_ADD_DUMMY_COUNTS \\
            --input_gene_results_dir "Rsem_gene_counts" \\
            --primary_keytype ${ meta.primary_keytype } \\
            --normalization 'ERCC-groupB' \\
            --normalized_counts_output_prefix "norm_counts_output/ERCC_" \\
            --dge_output_prefix "dge_output_ercc/ERCCnorm_" \\
            --annotation_file_path ${annotation_file} \\
            --extended_table_output_prefix "dge_output_ercc/"\\
            --extended_table_output_suffix "_ERCCnorm.csv" \\
            --verbose
    fi
    """

  script:
    """
    ./dge_annotation_R_scripts/dge_annotation_workflow.R \\
        --runsheet_path runsheet.csv \\
        --input_gene_results_dir "Rsem_gene_counts" \\
        --primary_keytype ${ meta.primary_keytype } \\
        --normalization 'default' \\
        --normalized_counts_output_prefix "norm_counts_output/" \\
        --dge_output_prefix "dge_output/" \\
        --annotation_file_path ${annotation_file} \\
        --extended_table_output_prefix "dge_output/"\\
        --extended_table_output_suffix ".csv" \\
        --verbose

    if ${ meta.has_ercc ? 'true' : 'false'}
    then
        ./dge_annotation_R_scripts/dge_annotation_workflow.R \\
            --runsheet_path runsheet.csv \\
            --input_gene_results_dir "Rsem_gene_counts" \\
            --primary_keytype ${ meta.primary_keytype } \\
            --normalization 'ERCC-groupB' \\
            --normalized_counts_output_prefix "norm_counts_output/ERCC_" \\
            --dge_output_prefix "dge_output_ercc/ERCCnorm_" \\
            --annotation_file_path ${annotation_file} \\
            --extended_table_output_prefix "dge_output_ercc/"\\
            --extended_table_output_suffix "_ERCCnorm.csv" \\
            --verbose
    fi
    """
}
