nextflow.enable.dsl=2
import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("MMddyyyy-HH_mm_ss")

c_bright_green = "\u001b[32;1m";
c_reset = "\033[0m";

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

workflow {
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

    meta_ch | DOWNLOAD_GENOME_ANNOTATIONS | set { genome_annotations }

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

    ALIGN_STAR.out | combine( BUILD_RSEM.out ) | set { aligned_ch }
    aligned_ch | COUNT_ALIGNED

    COUNT_ALIGNED.out | map { it[1] } | collect | set { rsem_ch }

    organism_ch = channel.fromPath( params.organismCSV )

    DGE_BY_DESEQ2( isa_ch, organism_ch, rsem_ch, meta_ch  )

    // VV processes
    ch_vv_log_00 = Channel.fromPath("nextflow_vv_log.tsv")
    VV_RAW_READS( raw_reads_ch | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes val(meta) from tuple
                  ch_vv_log_00 ) | set { ch_vv_log_01 }

    VV_RAW_READS_MULTIQC( RAW_MULTIQC.out.data,
                          ch_vv_log_01 ) | set { ch_vv_log_02 }

    VV_TRIMMED_READS( TRIMGALORE.out.reads | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes val(meta) from tuple
                      ch_vv_log_02 ) | set { ch_vv_log_03 }

    VV_TRIMMED_READS_MULTIQC( TRIMMED_MULTIQC.out.data,
                              ch_vv_log_03 ) | set { ch_vv_log_04 }

    VV_STAR_ALIGNMENTS( ALIGN_STAR.out | map{ it -> it[1..it.size()-1] } | collect, // map use here: removes val(meta) from tuple
                        ch_vv_log_04 ) | set { ch_vv_log_05 }

    VV_RSEM_COUNTS( COUNT_ALIGNED.out | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes val(meta) from tuple
                    ch_vv_log_05 ) | set { ch_vv_log_06 }

    VV_DESEQ2_ANALYSIS( DGE_BY_DESEQ2.out.dge | map{ it -> it[1..it.size()-1] } | collect, // map use here: removes val(meta) from tuple
                        ch_vv_log_06 ) | set { ch_vv_log_07 }

}

workflow.onComplete {
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Raw and Processed data location: ${ params.gldsAccession }"
    println "V&V logs location: ${ params.gldsAccession }-VV/VV_Log${c_reset}"
}
