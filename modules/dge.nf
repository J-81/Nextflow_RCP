/*
 * Different Gene Expression Analysis Processes
 */

/*
 * Script written by NASA, parameterized by dataset
 */

// TODO: REMOVE HARD CODED CONDA PATH

process DGE_BY_DESEQ2 {
  conda "${baseDir}/envs/RNAseq_Rtools.yml"
  publishDir "${ params.gldsAccession }/${meta.DESeq2_NormCount}", pattern: "norm_counts_output/*", saveAs: { "${file(it).getName()}" }
  publishDir "${ params.gldsAccession }/${meta.DESeq2_DGE}", pattern: "dge_output/*", saveAs: { "${file(it).getName()}" }
  if (meta.has_ercc) {
    publishDir "${ params.gldsAccession }/${meta.DESeq2_DGE}/ERCC_NormDGE", pattern: "dge_output_ercc/*", saveAs: { "${file(it).getName()}" }
  }

  input:
    path(Isa_zip)
    path(organisms_csv)
    path("Rsem_gene_counts/*")
    val(meta)
  output:
    tuple path("norm_counts_output/Normalized_Counts.csv"),
          path("norm_counts_output/SampleTable.csv"),
          path("norm_counts_output/ERCC_Normalized_Counts.csv"),
          path("norm_counts_output/Unnormalized_Counts.csv"), emit: norm_counts
    tuple path("dge_output/contrasts.csv"),
          path("dge_output/differential_expression.csv"),
          path("dge_output/visualization_output_table.csv"),
          path("dge_output/visualization_PCA_table.csv"), emit: dge
    path("dge_output_ercc"), emit: dge_ercc

    path("versions.txt"), emit: version
  script:
    def deseq2_script = meta.has_ercc ? "deseq2_normcounts_wERCC_DGE_vis_ISA.R" : "deseq2_normcounts_noERCC_DGE_vis_ISA.R"
    """
    # create output directories
    mkdir norm_counts_output
    mkdir dge_output
    mkdir dge_output_ercc

    # run the script with R
    ${deseq2_script} \
      ${ meta.organism_non_sci } \
      $Isa_zip \
      norm_counts_output \
      dge_output \
      dge_output_ercc
    """
}
