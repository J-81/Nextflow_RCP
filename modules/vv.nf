/* VV check processes
* These processes intentional deviate from Nextflow isolation to ensure QC reports
*   are based on files in publish directory and not work directories.
*/

// NOTE: first VV step also creates initial VV file that is shared across all vv steps
process VV_RAW_READS {
  // Log publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern:  "VV_log.tsv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process }.tsv" }
  // V&V'ed data publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '{00-RawData/Fastq,Metadata}',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_INPUT/Metadata/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_INPUT/00-RawData/Fastq/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_INPUT/00-RawData/FastQC_Reports/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_INPUT/00-RawData/FastQC_Reports/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC

  output:
    path("Metadata/*_runsheet.csv"), emit: VVed_runsheet
    path("00-RawData/Fastq"), emit: VVed_raw_reads
    path("00-RawData/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"), emit: VVed_raw_fastqc
    path("00-RawData/FastQC_Reports/raw_multiqc_report.zip"), emit: VVed_raw_multiqc_report
    path("VV_log.tsv"), optional: params.skipVV, emit: log

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !params.skipVV} ; then
      raw_reads_VV.py  --root-path . --accession ${ params.gldsAccession }
    fi
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
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !params.skipVV} ; then
      trimmed_reads_VV.py --root-path . --accession ${ params.gldsAccession }
    fi
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