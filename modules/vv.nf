/* VV check processes
*/

// NOTE: first VV step also creates inital VV file that is shared across all vv steps
process VV_RAW_READS {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    val(samples)
    path(raw_reads)
    path(vv_config)

  output:
    path("VV_log_raw.txt")

  script:
    """
    raw_reads_VV.py --config ${ vv_config } \
                    --samples ${ samples.join(' ') } \
                    --input ${ raw_reads } \
                    --output VV_log_raw.txt
    """
}

process VV_RAW_READS_MULTIQC {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    val(samples)
    path(multiqcDataDir)
    path(vv_config)
    path(vv_output)

  output:
    path("VV_log_raw_read_multiqc.txt")

  script:
    """
    raw_reads_multiqc_VV.py --config ${ vv_config } \
                            --samples ${ samples.join(' ') } \
                            --input ${ multiqcDataDir } \
                            --output ${ vv_output }
    cp -L ${ vv_output } VV_log_raw_read_multiqc.txt
    """
}


process VV_TRIMMED_READS {
  stageInMode "copy"
  // publishDir "${params.publishDirPath}/VV/${params.timestamp}",

  input:
    val(samples)
    path(trimmed_reads)
    path(vv_config)
    path(vv_output)

  output:
    path("VV_log_trimmed.txt")

  script:
    """
    trimmed_reads_VV.py --config ${ vv_config } \
                    --samples ${ samples.join(' ') } \
                    --input ${ trimmed_reads } \
                    --output ${ vv_output }
    cp -L ${ vv_output } VV_log_trimmed.txt
    """
}

process VV_TRIMMED_READS_MULTIQC {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    val(samples)
    path(multiqcDataDir)
    path(vv_config)
    path(vv_output)

  output:
    path("VV_log_trimmed_read_multiqc.txt")

  script:
    """
    trimmed_reads_multiqc_VV.py --config ${ vv_config } \
                                --samples ${ samples.join(' ') } \
                                --input ${ multiqcDataDir } \
                                --output ${ vv_output }
    cp -L ${ vv_output } VV_log_trimmed_read_multiqc.txt
    """
}

process VV_STAR_ALIGNMENTS {
  stageInMode "copy"
  publishDir "${params.publishDirPath}/VV/${params.timestamp}",
              mode: 'copy', saveAs: { "VV_RESULTS.txt" }


  input:
    val(samples)
    path(genomeMapping)
    path(transcriptomeMapping)
    path(logs)
    path(vv_config)
    path(vv_output)

  output:
    path("VV_log_star_alignments.txt")

  script:
    """
    star_alignments_VV.py --config ${ vv_config } \
                          --samples ${ samples.join(' ') } \
                          --g ${ genomeMapping } \
                          --t ${ transcriptomeMapping } \
                          --l ${ logs } \
                          --output ${ vv_output }
    cp -L ${ vv_output } VV_log_star_alignments.txt
    """
}
