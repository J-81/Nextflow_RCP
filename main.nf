nextflow.enable.dsl=2
import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("MMddyyyy-HH_mm_ss")

include { DOWNLOAD_RAW_READS;
          DOWNLOAD_GENOME_ANNOTATIONS;
          DOWNLOAD_ERCC;
          DOWNLOAD_ISA } from './modules/download.nf'
include { FASTQC as RAW_FASTQC } from './modules/quality.nf' addParams(fastQCLabel: 'raw')
include { FASTQC as TRIM_FASTQC } from './modules/quality.nf' addParams(fastQCLabel: 'trimmed')
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'raw')
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'trimmed')
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_STAR;
          ALIGN_STAR;
          BUILD_RSEM;
          COUNT_ALIGNED;
          SUBSAMPLE_GENOME;
          CONCAT_ERCC } from './modules/genome.nf'
include { DGE_BY_DESEQ2 } from './modules/dge.nf'
include { SAMPLES_FROM_ISA } from './modules/isa.nf'
include { VV_RAW_READS;
          VV_TRIMMED_READS;
          VV_RAW_READS_MULTIQC;
          VV_TRIMMED_READS_MULTIQC;
          VV_STAR_ALIGNMENTS } from './modules/vv.nf' addParams(timestamp: sdf.format(date))

/*
 * Starting point, includes downloads data from GeneLab
 */

workflow GET_DATA {
  take: samples_ch
  main:
    samples_ch |  DOWNLOAD_RAW_READS
  emit:
    DOWNLOAD_RAW_READS.out.raw_reads
}

workflow {
	main:
    DOWNLOAD_ISA | set{ isa_ch }
    SAMPLES_FROM_ISA( isa_ch ) | splitText { it.replaceAll("\\s","") } | set{ samples_ch }

    if ( params.raw_reads ) {
      raw_reads_ch = Channel.from( params.raw_reads )
                            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
    } else {
      GET_DATA( samples_ch ) | set { raw_reads_ch }
    }

    if ( params.vv_config_file ) {
      raw_reads_ch | multiMap { it ->
                      samples: it[0]
                      input_files: it[1]
                      }
                   | set{ raw_reads_vv_input }

      VV_RAW_READS( raw_reads_vv_input.samples | collect,
                    raw_reads_vv_input.input_files | collect,
                    params.vv_config_file ) | set { vv_output_ch }
    }

    raw_reads_ch | RAW_FASTQC

    RAW_FASTQC.out | map { it -> [ it[1], it[2] ] } \
                   | flatten \
                   | unique \
                   | collect \
                   | RAW_MULTIQC

    if ( params.vv_config_file ) {
      VV_RAW_READS_MULTIQC( samples_ch | collect,
                            RAW_MULTIQC.out.data,
                            params.vv_config_file,
                            vv_output_ch ) | set { vv_output_ch }
    }

    raw_reads_ch | map{ it -> [ it[0], it[1][0], it[1][1] ] } | TRIMGALORE

    if ( params.vv_config_file ) {
      TRIMGALORE.out.reads | multiMap { it ->
                              samples: it[0]
                              input_files: it[1]
                              }
                           | set{ trimmed_reads_vv_input }

      VV_TRIMMED_READS( trimmed_reads_vv_input.samples | collect,
                        trimmed_reads_vv_input.input_files | collect,
                        params.vv_config_file,
                        vv_output_ch ) | set { vv_output_ch }
    }

    TRIMGALORE.out.reads | map{ it -> [ it[0], [ it[1], it[2] ] ]} | TRIM_FASTQC

    TRIM_FASTQC.out | map { it -> [ it[1], it[2] ] } \
                    | flatten \
                    | unique \
                    | collect \
                    | TRIM_MULTIQC

    if ( params.vv_config_file ) {
      VV_TRIMMED_READS_MULTIQC( samples_ch | collect,
                                TRIM_MULTIQC.out.data,
                                params.vv_config_file,
                                vv_output_ch ) | set { vv_output_ch }
    }

    if ( params.genomeFasta && params.genomeGTF ) {
      genome_annotations = channel.fromPath([ params.genomeFasta, params.genomeGTF ])
                                  .toList()
    } else {
      DOWNLOAD_GENOME_ANNOTATIONS | set { genome_annotations }
    }

    /*
    SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
    */
    if ( params.genomeSubsample ) {
      genome_annotations | SUBSAMPLE_GENOME | set { genome_annotations }
    }

    /*
    ERCC STEP : ADD ERCC Fasta and GTF to genome files
    */
    if ( params.ERCC ) {
      DOWNLOAD_ERCC | set { ercc_annotations }
      CONCAT_ERCC( genome_annotations, ercc_annotations ) | set { genome_annotations }
    }

    genome_annotations | BUILD_STAR

    TRIMGALORE.out.reads | combine( BUILD_STAR.out ) | ALIGN_STAR

    if ( params.vv_config_file ) {
      VV_STAR_ALIGNMENTS( samples_ch | collect,
                          ALIGN_STAR.out.genomeMapping | map{ it -> it[1] } | collect,
                          ALIGN_STAR.out.transcriptomeMapping | map{ it -> it[1] } | collect,
                          ALIGN_STAR.out.logs| map{ it -> it[1..it.size()-1] } | collect,
                          params.vv_config_file,
                          vv_output_ch ) | set { vv_output_ch }
    }

    genome_annotations | BUILD_RSEM

    ALIGN_STAR.out.transcriptomeMapping | combine( BUILD_RSEM.out ) | set { aligned_ch }
    aligned_ch | COUNT_ALIGNED

    COUNT_ALIGNED.out.countsPerGene | map { it[1] } | collect | toList | set { rsem_ch }

    organism_ch = channel.fromPath( params.organismCSV )
    external_ch = isa_ch.combine( organism_ch )

    external_ch | combine( rsem_ch ) | set { for_dge_ch }
    for_dge_ch | DGE_BY_DESEQ2

}
