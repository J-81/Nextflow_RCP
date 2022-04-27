/* VV check processes
* These processes intentional deviate from Nextflow isolation to ensure QC reports
*   are based on files in publish directory and not work directories.
*/

// NOTE: first VV step also creates initial VV file that is shared across all vv steps
process VV_RAW_READS {
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_verbose_through_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path(vv_log)

  output:
    path("VV_out.tsv")

  script:
    """
    raw_reads_VV.py  --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    mv VV_log_verbose.tsv VV_out.tsv
    """
}

process VV_TRIMMED_READS {
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_verbose_through_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    trimmed_reads_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    tail -n +2 VV_log_verbose.tsv >> VV_in.tsv
    mv VV_in.tsv VV_out.tsv
    """
}


process VV_STAR_ALIGNMENTS {
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_verbose_through_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    star_alignments_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    tail -n +2 VV_log_verbose.tsv >> VV_in.tsv
    mv VV_in.tsv VV_out.tsv
    """
}

process VV_RSEQC {
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_verbose_through_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    rseqc_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    tail -n +2 VV_log_verbose.tsv >> VV_in.tsv
    mv VV_in.tsv VV_out.tsv
    """
}


process VV_RSEM_COUNTS {
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode,
    saveAs: { "VV_log_verbose_through_${ task.process }.tsv" }

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    rsem_counts_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    tail -n +2 VV_log_verbose.tsv >> VV_in.tsv
    mv VV_in.tsv VV_out.tsv
    """
}

process VV_DESEQ2_ANALYSIS {
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.RootDirForVV }/VV_Logs",
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    tuple path("VV_log_final.tsv"), path("VV_log_final_only_issues.tsv")

  stub:
    // SET MAX FLAG CODE TO ONLY HALT ON DEVELOPER LEVEL FLAGS
    """
    deseq2_script_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession } --max-flag-code 90
    tail -n +2 VV_log_verbose.tsv >> VV_in.tsv
    mv VV_in.tsv VV_log_final.tsv

    # filtered log
    filer_to_only_issues.py
    """

  script:
    """
    deseq2_script_VV.py --root-path ${ params.RootDirForVV } --accession ${ params.gldsAccession }
    tail -n +2 VV_log_verbose.tsv >> VV_in.tsv
    mv VV_in.tsv VV_log_final.tsv

    # filtered log
    filer_to_only_issues.py
    """
}
