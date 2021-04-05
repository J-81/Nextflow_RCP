workflow VV {
  take:
    meta_ch = meta_ch
    // 00-RawData
    raw_reads_fastq                          = raw_reads_ch
    raw_reads_fastqc                         = RAW_FASTQC.out
    raw_reads_multiqc                        = RAW_MULTIQC.out
    // 01-TG_Preproc
    trimmed_reads_fastq                      = TRIMGALORE.out.reads
    trimmed_reads_fastqc                     = TRIMMED_FASTQC.out
    trimmed_reads_multiqc                    = TRIMMED_MULTIQC.out
    // 02-STAR_Alignment
    star_alignments_genomeMapping            = ALIGN_STAR.out.genomeMapping
    star_alignments_transcriptomeMapping     = ALIGN_STAR.out.transcriptomeMapping
    star_alignments_logs                     = ALIGN_STAR.out.logs
    // 03-RSEM_Counts
    rsem_counts_countsPerGene                = COUNT_ALIGNED.out.countsPerGene
    rsem_counts_countsPerIsoform             = COUNT_ALIGNED.out.countsPerIsoform
    rsem_counts_stats                        = COUNT_ALIGNED.out.stats
    // 04-DESeq2_NormCounts & 05-DESeq2_DGE
    deseq2_normCounts                        = DGE_BY_DESEQ2.out.norm_counts,
    deseq2_dge                               = DGE_BY_DESEQ2.out.dge

    main:
  // Validation and Verification Block
  // - Nextflow Note: Even though this is defined at the end, Nextflow will
  //     execute VV tasks once the required input files are ready. i.e. VV is
  //     performed in parallel with task processing.
    ch_vv_log = Channel.fromPath("nextflow_vv_log.tsv")
    VV_RAW_READS( meta_ch,
                  raw_reads_fastq | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes value sample from tuple
                  params.vv_config_file ) | mix(vv_output_ch) | set{ vv_output_ch }
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
