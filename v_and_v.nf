workflow VV {
  take:
    raw_reads_ch
  main:
      // Validation and Verification Block
      // - Nextflow Note: Even though this is defined at the end, Nextflow will
      //     execute VV tasks once the required input files are ready. i.e. VV is
      //     performed in parallel with task processing.
        raw_reads_ch | view
        ch_vv_log = Channel.fromPath("nextflow_vv_log.tsv")
        VV_RAW_READS( raw_reads_ch,
                      raw_reads_ch | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes value sample from tuple
                      ch_vv_log ) | mix(vv_output_ch) | set{ vv_output_ch }
      /*
      VV_RAW_READS_MULTIQC( samples_ch | collectFile(name: "samples.txt", newLine: true),
                            RAW_MULTIQC.out.data,
                            params.vv_config_file) | mix(vv_output_ch) | set{ vv_output_ch }

      VV_TRIMMED_READS( samples_ch | collectFile(name: "samples.txt", newLine: true),
                        TRIMGALORE.out.reads | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes value sample from tuple
                        params.vv_config_file)  | mix(vv_output_ch) | set{ vv_output_ch }

      VV_TRIMMED_READS_MULTIQC( samples_ch | collectFile(name: "samples.txt", newLine: true),
                                TRIM_MULTIQC.out.data,
                                params.vv_config_file)  | mix(vv_output_ch) | set{ vv_output_ch }

      VV_STAR_ALIGNMENTS( samples_ch | collectFile(name: "samples.txt", newLine: true),
                          ALIGN_STAR.out.genomeMapping | map{ it -> it[1] } | collect,
                          ALIGN_STAR.out.transcriptomeMapping | map{ it -> it[1] } | collect,
                          ALIGN_STAR.out.logs| map{ it -> it[1..it.size()-1] } | collect,
                          params.vv_config_file) | mix(vv_output_ch) | set{ vv_output_ch }

      VV_RSEM_COUNTS( samples_ch | collectFile(name: "samples.txt", newLine: true),
                      COUNT_ALIGNED.out.countsPerGene | map{ it -> it[1] } | collect,
                      COUNT_ALIGNED.out.countsPerIsoform | map{ it -> it[1] } | collect,
                      COUNT_ALIGNED.out.stats | map{ it -> it[1] } | collect,
                      params.vv_config_file)  | mix(vv_output_ch) | set{ vv_output_ch }

      VV_DESEQ2_ANALYSIS( samples_ch | collectFile(name: "samples.txt", newLine: true),
                          DGE_BY_DESEQ2.out.norm_counts,
                          DGE_BY_DESEQ2.out.dge,
                          params.vv_config_file) | mix(vv_output_ch) | set{ vv_output_ch }
      vv_output_ch | map{ it.text } | collectFile(name: "${sdf.format(date)}_Final_VV.tsv", storeDir: workflow.launchDir) | view
      */
}
