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
  publishDir "${ params.gldsAccession }/${meta.DESeq2_DGE}/ERCC_NormDGE", pattern: "dge_output_ercc/*", saveAs: { "${file(it).getName()}" }

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
    tuple path("dge_output_ercc/ERCCnorm_contrasts.csv"),
          path("dge_output_ercc/ERCCnorm_differential_expression.csv"),
          path("dge_output_ercc/visualization_output_table_ERCCnorm.csv"),
          path("dge_output_ercc/visualization_PCA_table_ERCCnorm.csv"), optional: !params.ERCC, emit: dge_ercc

    path("versions.txt"), emit: version
  script:
    """
    # create output directories
    mkdir norm_counts_output
    mkdir dge_output
    mkdir dge_output_ercc

    # run the script with R
    deseq2_normcounts_wERCC_DGE_vis_ISA.R \
      ${ meta.organism_non_sci } \
      $Isa_zip \
      norm_counts_output \
      dge_output \
      dge_output_ercc

    Rscript -e "library(tximport); write(x=as.character(packageVersion('tximport')), file='version.txt')"
    Rscript -e "library(DESeq2); write(x=as.character(packageVersion('DESeq2')), file='version.txt')"
    Rscript -e "library(tidyverse); write(x=as.character(packageVersion('tidyverse')), file='version.txt')"
    Rscript -e "library(Risa); write(x=as.character(packageVersion('Risa')), file='version.txt')"
    Rscript -e "library(STRINGdb); write(x=as.character(packageVersion('STRINGdb')), file='version.txt')"
    Rscript -e "library(PANTHER.db); write(x=as.character(packageVersion('PANTHER.db')), file='version.txt')"
    # TODO: add version for annotations db (organism dependent) Rscript -e "library(PANTHER.db); write(x=as.character(packageVersion('PANTHER.db')), file='version.txt')"
    """
}
