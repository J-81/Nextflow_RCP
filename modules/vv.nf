/* VV check processes
*/

// NOTE: first VV step also creates inital VV file that is shared across all vv steps
process VV_RAW_READS {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    path(samples)
    path(raw_reads), stageAs: "rawReads/*"
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    raw_reads_VV.py --config ${ vv_config } \
                    --samples ${ samples } \
                    --input rawReads \
                    --output VV_out.tsv
    """
}

process VV_RAW_READS_MULTIQC {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    val(samples)
    path(multiqcDataDir)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    raw_reads_multiqc_VV.py --config ${ vv_config } \
                            --samples ${ samples.join(' ') } \
                            --input ${ multiqcDataDir } \
                            --output VV_out.tsv
    """
}


process VV_TRIMMED_READS {
  stageInMode "copy"
  // publishDir "${params.publishDirPath}/VV/${params.timestamp}",

  input:
    val(samples)
    path(trimmed_reads)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    trimmed_reads_VV.py --config ${ vv_config } \
                    --samples ${ samples.join(' ') } \
                    --input ${ trimmed_reads } \
                    --output VV_out.tsv
    """
}

process VV_TRIMMED_READS_MULTIQC {
  stageInMode "copy"
  //publishDir "${params.publishDirPath}/VV/${params.timestamp}", mode: 'copy'

  input:
    val(samples)
    path(multiqcDataDir)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    trimmed_reads_multiqc_VV.py --config ${ vv_config } \
                                --samples ${ samples.join(' ') } \
                                --input ${ multiqcDataDir } \
                                --output VV_out.tsv
    """
}

process VV_STAR_ALIGNMENTS {
  stageInMode "copy"
  // publishDir "${params.publishDirPath}/VV/${params.timestamp}",
  //             mode: 'copy', saveAs: { "VV_RESULTS.txt" }


  input:
    val(samples)
    path(genomeMapping)
    path(transcriptomeMapping)
    path(logs)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    star_alignments_VV.py --config ${ vv_config } \
                          --samples ${ samples.join(' ') } \
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
    val(samples)
    path(geneCounts)
    path(transcriptCounts)
    path(stats)
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    rsem_counts_VV.py --config ${ vv_config } \
                      --samples ${ samples.join(' ') } \
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
    val(samples)
    path(norm_countsDir), stageAs: 'counts/*'
    path(dgeDir), stageAs: 'dge/*'
    path(vv_config)

  output:
    path("VV_out.tsv")

  script:
    """
    deseq2_script_VV.py --config ${ vv_config } \
                        --samples ${ samples.join(' ') } \
                        --normDir counts \
                        --dgeDir dge \
                        --output VV_out.tsv
    """
}
