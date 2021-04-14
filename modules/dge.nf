/*
 * Different Gene Expression Analysis Processes
 */

/*
 * Script written by NASA, parameterized by dataset
 */

// TODO: REMOVE HARD CODED CONDA PATH

process DGE_BY_DESEQ2 {
  //conda "${baseDir}/envs/RNAseq_Rtools.yml"
  publishDir "${ params.gldsAccession }" }
  // publishDir "${ params.gldsAccession }/${meta.DESeq2_DGE}", pattern: "dge_output/*", saveAs: { "${file(it).getName()}" }
  // publishDir "${ params.gldsAccession }/${meta.DESeq2_DGE}/ERCC_NormDGE", pattern: "dge_output_ercc/*", saveAs: { "${file(it).getName()}" }

  input:
    path(Isa_zip)
    path(organisms_csv)
    path("Rsem_gene_counts/*")
    val(meta)
  output:
    path("modified_dge.txt")
  script:
    def deseq2_script = meta.has_ercc ? "deseq2_normcounts_wERCC_DGE_vis_ISA.R" : "deseq2_normcounts_noERCC_DGE_vis_ISA.R"
    """
    echo Now I just announce that the presenter would be happy to take questions > modified_dge.txt
    """
}
