/* VV check processes
* These processes intentional deviate from Nextflow isolation to ensure QC reports
*   are based on files in publish directory and not work directories.
*/

// NOTE: first VV step also creates initial VV file that is shared across all vv steps
process VV_RAW_READS {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path(vv_log)

  output:
    path("VV_out.tsv")

  script:
    """
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    # ensure no existing VV_out.tsv file
    rm -rf VV_Log

    raw_reads_VV.py  --runsheet-path Metadata/*runsheet.csv \
                     --output VV_Log_Up_To_Halt.tsv \
                     --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
    """
}

process VV_RAW_READS_MULTIQC {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    raw_reads_multiqc_VV.py --runsheet-path Metadata/*runsheet.csv \
                            --output VV_Log_Up_To_Halt.tsv \
                            --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
    """
}


process VV_TRIMMED_READS {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    trimmed_reads_VV.py --runsheet-path Metadata/*runsheet.csv \
                        --output VV_Log_Up_To_Halt.tsv \
                        --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
    """
}

process VV_TRIMMED_READS_MULTIQC {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
  """
  # copy to processed data directory
  cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv
  # cd into processed data directory
  cd ${ params.RootDirForVV }/${ params.gldsAccession }
  trimmed_reads_multiqc_VV.py --runsheet-path Metadata/*runsheet.csv \
                              --output VV_Log_Up_To_Halt.tsv \
                              --halt-severity 90
  # move back to work dir and mv tsv into work dir
  cd -
  mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
  """
}

process VV_STAR_ALIGNMENTS {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    star_alignments_VV.py --runsheet-path Metadata/*runsheet.csv \
                          --output VV_Log_Up_To_Halt.tsv \
                          --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
    """
}

process VV_RSEQC {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    rseqc_VV.py --runsheet-path Metadata/*runsheet.csv \
                          --output VV_Log_Up_To_Halt.tsv \
                          --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
    """
}


process VV_RSEM_COUNTS {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    rsem_counts_VV.py --runsheet-path Metadata/*runsheet.csv \
                      --output VV_Log_Up_To_Halt.tsv \
                      --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log_Up_To_Halt.tsv VV_out.tsv
    """
}

process VV_DESEQ2_ANALYSIS {
  tag "Dataset: ${ params.gldsAccession }"

  label 'VV'

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_Log")

  stub:
    """
    # copy to processed data directory
    mkdir -p ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log/VV_FULL_OUT.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    deseq2_script_VV.py --runsheet-path Metadata/*runsheet.csv \
                        --output VV_Log/VV_FULL_OUT.tsv \
                        --halt-severity 91 # required as the stub deseq2 script results in counts that differ from the rsem, a halting error by default

    # move back to work dir and mv tsv into work dir
    cd -
    # mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log
    touch VV_Log # signals end of VV
    """

  script:
    """
    # copy to processed data directory
    mkdir -p ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log
    cp -L VV_in.tsv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log/VV_FULL_OUT.tsv
    # cd into processed data directory
    cd ${ params.RootDirForVV }/${ params.gldsAccession }
    deseq2_script_VV.py --runsheet-path Metadata/*runsheet.csv \
                        --output VV_Log/VV_FULL_OUT.tsv \
                        --halt-severity 90

    # move back to work dir and mv tsv into work dir
    cd -
    # mv ${ params.RootDirForVV }/${ params.gldsAccession }/VV_Log
    touch VV_Log # signals end of VV
    """
}
