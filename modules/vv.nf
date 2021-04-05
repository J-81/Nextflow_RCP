/* VV check processes
*/

// NOTE: first VV step also creates inital VV file that is shared across all vv steps
process VV_RAW_READS {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    val(meta)
    path(runsheet)
    path("${ meta.raw_read_root_dir }/*")
    path(vv_log)

  output:
    path("VV_out.tsv")

  script:
    """
    # Not needed for first VV process
    # cp -L ${vv_log} appendTo.tsv
    raw_reads_VV.py  --runsheet-path ${ runsheet } \
                     --output appendTo.tsv \
                     --halt-severity 90
    mv appendTo.tsv VV_out.tsv
    """
}

process VV_RAW_READS_MULTIQC {
  conda "${baseDir}/envs/VV.yml"
  publishDir "VV"

  input:
    val(meta)
    path(runsheet)
    // resulted in nested staging:
    // path("${ meta.raw_read_multiqc }/*")
    path("00-RawData/FastQC_Reports/raw_multiqc_report/*")
    path(multiqcHtmlPath)
    path("VV_in.tsv")

  output:
    path("VV_out.tsv")

  script:
    """
    cp -L VV_in.tsv appendTo.tsv
    raw_reads_multiqc_VV.py --runsheet-path ${ runsheet } \
                            --output appendTo.tsv \
                            --halt-severity 90
    cp appendTo.tsv VV_out.tsv
    """
}


process VV_TRIMMED_READS {
  stageInMode "copy"
  // publishDir "${params.publishDirPath}/VV/${params.timestamp}",

  input:
    path(samples)
    path(trimmed_reads)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    trimmed_reads_VV.py --config ${ vv_config } \
                    --samples ${ samples } \
                    --input ${ trimmed_reads } \
                    --output VV_out.tsv
    """
}

process VV_TRIMMED_READS_MULTIQC {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    path(samples)
    path(multiqcDataDir)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    trimmed_reads_multiqc_VV.py --config ${ vv_config } \
                                --samples ${ samples } \
                                --input ${ multiqcDataDir } \
                                --output VV_out.tsv
    """
}

process VV_STAR_ALIGNMENTS {
  stageInMode "copy"
  // publishDir "${params.publishDirPath}/VV/${params.timestamp}",
  //             mode: 'copy', saveAs: { "VV_RESULTS.txt" }


  input:
    path(samples)
    path(genomeMapping)
    path(transcriptomeMapping)
    path(logs)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    star_alignments_VV.py --config ${ vv_config } \
                          --samples ${ samples } \
                          --g ${ genomeMapping } \
                          --t ${ transcriptomeMapping } \
                          --l ${ logs } \
                          --output VV_out.tsv
    """
}

process VV_RSEM_COUNTS {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}",
  //            mode: 'copy', saveAs: { "VV_RESULTS.txt" }


  input:
    path(samples)
    path(geneCounts)
    path(transcriptCounts)
    path(stats)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    rsem_counts_VV.py --config ${ vv_config } \
                      --samples ${ samples } \
                      --g ${ geneCounts } \
                      --t ${ transcriptCounts } \
                      --stats ${ stats } \
                      --output VV_out.tsv
    """
}

process VV_DESEQ2_ANALYSIS {
  stageInMode "copy"
  publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    path(samples)
    path(norm_countsDir), stageAs: 'counts/*'
    path(dgeDir), stageAs: 'dge/*'
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    deseq2_script_VV.py --config ${ vv_config } \
                        --samples ${ samples } \
                        --normDir counts \
                        --dgeDir dge \
                        --output VV_out.tsv
    """
}
