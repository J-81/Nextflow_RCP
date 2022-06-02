/* VV check processes
* These processes intentional deviate from Nextflow isolation to ensure QC reports
*   are based on files in publish directory and not work directories.
*/

// NOTE: first VV step also creates initial VV file that is shared across all vv steps
process VV_RAW_READS {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path(vv_log)

  output:
    path("VV_log.tsv")

  script:
    """
    raw_reads_VV.py  --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    """
}

process VV_TRIMMED_READS {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_log.tsv")

  script:
    """
    trimmed_reads_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    """
}


process VV_STAR_ALIGNMENTS {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_log.tsv")

  script:
    """
    star_alignments_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    """
}

process VV_RSEQC {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_log.tsv")

  script:
    """
    rseqc_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    """
}


process VV_RSEM_COUNTS {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_log.tsv")

  script:
    """
    rsem_counts_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    """
}

process VV_DESEQ2_ANALYSIS {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_log.tsv")

  stub:
    // SET MAX FLAG CODE TO ONLY HALT ON DEVELOPER LEVEL FLAGS
    """
    deseq2_script_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession } --max-flag-code 90
    """

  script:
    """
    deseq2_script_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    """
}

process VV_CONCAT_FILTER {
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_in.tsv")

  output:
    tuple path("VV_log_final.tsv"), path("VV_log_final_only_issues.tsv")

  script:
    """
    concat_logs.py
    filter_to_only_issues.py
    """
}