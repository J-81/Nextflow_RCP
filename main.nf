nextflow.enable.dsl=2
import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("MMddyyyy-HH_mm_ss")

include { DOWNLOAD_RAW_READS;
          DOWNLOAD_GENOME_ANNOTATIONS;
          DOWNLOAD_ERCC;
          DOWNLOAD_ISA } from './modules/download.nf'
include { RAW_FASTQC
          TRIMMED_FASTQC
          RAW_MULTIQC
          TRIMMED_MULTIQC } from './modules/quality.nf'
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_STAR;
          ALIGN_STAR;
          BUILD_RSEM;
          COUNT_ALIGNED;
          SUBSAMPLE_GENOME;
          CONCAT_ERCC } from './modules/genome.nf'
include { DGE_BY_DESEQ2 } from './modules/dge.nf'
include { PARSE_ISA } from './modules/isa.nf'
include { VV_RAW_READS;
          VV_TRIMMED_READS;
          VV_RAW_READS_MULTIQC;
          VV_TRIMMED_READS_MULTIQC;
          VV_STAR_ALIGNMENTS;
          VV_RSEM_COUNTS;
          VV_DESEQ2_ANALYSIS } from './modules/vv.nf' addParams(timestamp: sdf.format(date))

include { staging as STAGING } from './stage_analysis.nf'
println "PARAMS: $params"

workflow RNASEQ {
	main:
    STAGING()
    STAGING.out.runsheet
    STAGING.out.raw_reads | set { raw_reads_ch }
    // meta only for dataset specific processes that don't use samples
    // e.g. downloading correct reference genome base on organism
    STAGING.out.raw_reads | take(1) | map{it -> it[0]} | view {"META: $it"} | set { meta_ch }
    STAGING.out.isa | set { isa_ch }

    raw_reads_ch | RAW_FASTQC //| view {"POST_FASTQC: $it"}

    RAW_FASTQC.out | map { it -> [ it[1], it[2] ] }
                   //| view {"POST_MAP: $it"}
                   | flatten
                   | unique
                   | collect
                   //| view {"PRE_RAW_MULTIQC: $it"}
                   | RAW_MULTIQC

    raw_reads_ch |  TRIMGALORE

    TRIMGALORE.out.reads | TRIMMED_FASTQC

    TRIMMED_FASTQC.out | map { it -> [ it[1], it[2] ] } \
                       | flatten \
                       | unique \
                       | collect \
                       | TRIMMED_MULTIQC

    if ( params.genomeFasta && params.genomeGTF ) {
      genome_annotations = channel.fromPath([ params.genomeFasta, params.genomeGTF ])
                                  .toList()
    } else {
      meta_ch | DOWNLOAD_GENOME_ANNOTATIONS | set { genome_annotations }
    }

    // SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
    if ( params.genomeSubsample ) {
      SUBSAMPLE_GENOME( genome_annotations, meta_ch ) | set { genome_annotations }
    }

    // ERCC STEP : ADD ERCC Fasta and GTF to genome files
    if ( params.ERCC ) {
      DOWNLOAD_ERCC | set { ercc_annotations }
      CONCAT_ERCC( genome_annotations, ercc_annotations, meta_ch ) | set { genome_annotations }
    }

    BUILD_STAR( genome_annotations, meta_ch)

    TRIMGALORE.out.reads | combine( BUILD_STAR.out ) | ALIGN_STAR

    BUILD_RSEM( genome_annotations, meta_ch)

    ALIGN_STAR.out.transcriptomeMapping | combine( BUILD_RSEM.out ) | set { aligned_ch }
    aligned_ch | COUNT_ALIGNED

    COUNT_ALIGNED.out.countsPerGene | map { it[1] } | collect | set { rsem_ch }

    organism_ch = channel.fromPath( params.organismCSV )

    DGE_BY_DESEQ2( isa_ch, organism_ch, rsem_ch, meta_ch  )

    /*
    // Validation and Verification Block
    // - Nextflow Note: Even though this is defined at the end, Nextflow will
    //     execute VV tasks once the required input files are ready. i.e. VV is
    //     performed in parallel with task processing.
    if ( params.vv_config_file ) {
      VV_RAW_READS( samples_ch | collectFile(name: "samples.txt", newLine: true),
                    raw_reads_ch | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes value sample from tuple
                    params.vv_config_file ) | mix(vv_output_ch) | set{ vv_output_ch }

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

    }
    */
}
