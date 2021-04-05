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
  emit:
      meta_ch                                  = meta_ch
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
      */
}
