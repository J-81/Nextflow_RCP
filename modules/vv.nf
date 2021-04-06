/* VV check processes
*/

// NOTE: first VV step also creates inital VV file that is shared across all vv steps
process VV_RAW_READS {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path(vv_log)

  output:
    path("VV_out.tsv")

  script:
    """
    cd ${workflow.launchDir}/${ params.gldsAccession }

    raw_reads_VV.py  --runsheet-path Metadata/*runsheet.csv \
                     --output appendTo.tsv \
                     --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
    """
}

process VV_RAW_READS_MULTIQC {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv
    # cd into processed data directory
    cd ${workflow.launchDir}/${ params.gldsAccession }
    raw_reads_multiqc_VV.py --runsheet-path Metadata/*runsheet.csv \
                            --output appendTo.tsv \
                            --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
    """
}


process VV_TRIMMED_READS {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv
    # cd into processed data directory
    cd ${workflow.launchDir}/${ params.gldsAccession }
    trimmed_reads_VV.py --runsheet-path Metadata/*runsheet.csv \
                        --output appendTo.tsv \
                        --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
    """
}

process VV_TRIMMED_READS_MULTIQC {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
  """
  # copy to processed data directory
  cp -L VV_in.tsv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv
  # cd into processed data directory
  cd ${workflow.launchDir}/${ params.gldsAccession }
  trimmed_reads_multiqc_VV.py --runsheet-path Metadata/*runsheet.csv \
                              --output appendTo.tsv \
                              --halt-severity 90
  # move back to work dir and mv tsv into work dir
  cd -
  mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
  """
}

process VV_STAR_ALIGNMENTS {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv
    # cd into processed data directory
    cd ${workflow.launchDir}/${ params.gldsAccession }
    star_alignments_VV.py --runsheet-path Metadata/*runsheet.csv \
                          --output appendTo.tsv \
                          --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
    """
}

process VV_RSEM_COUNTS {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv
    # cd into processed data directory
    cd ${workflow.launchDir}/${ params.gldsAccession }
    rsem_counts_VV.py --runsheet-path Metadata/*runsheet.csv \
                      --output appendTo.tsv \
                      --halt-severity 90
    # move back to work dir and mv tsv into work dir
    cd -
    mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
    """
}

process VV_DESEQ2_ANALYSIS {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    path("NULL") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    # copy to processed data directory
    cp -L VV_in.tsv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv
    # cd into processed data directory
    cd ${workflow.launchDir}/${ params.gldsAccession }
    deseq2_script_VV.py --runsheet-path Metadata/*runsheet.csv \
                        --output appendTo.tsv \
                        --halt-severity 90
    # remove temporary log used by null flagger
    rm tmp_remove.tsv


    # move back to work dir and mv tsv into work dir
    cd -
    mv ${workflow.launchDir}/${ params.gldsAccession }/appendTo.tsv VV_out.tsv
    """
}
