nextflow.enable.dsl=2
import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("MMddyyyy-HH_mm_ss")
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { DOWNLOAD_RAW_READS;
          DOWNLOAD_GENOME_ANNOTATIONS;
          DOWNLOAD_ERCC;
          DOWNLOAD_ISA } from './modules/download.nf'
include { RAW_FASTQC
          TRIMMED_FASTQC
          RAW_MULTIQC
          TRIMMED_MULTIQC
          ALIGN_MULTIQC } from './modules/quality.nf'
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

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ RNASeq Concensus Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("usage: nextflow run J-81/Nextflow_RCP -r help_menu_update --gldsAccession GLDS-000 --ensemblVersion 96  [--stageLocal] [--limitSamplesTo n] [--truncateTo n] [--genomeSubsample n]")
  println()
  println("required arguments:")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession number to stage raw reads for the RNASeq Concensus Pipeline")
  println("  --ensemblVersion n    ensembl Version to use for the reference genome. Default: 96")
  println("optional arguments:")
  println("  --help                show this help message and exit")
  println("  --skipVV              skip automated V&V checks")
  println("  --limitSamplesTo n    limit the number of samples staged to a number.")
  println("  --genomeSubsample n   subsamples genome fasta and gtf files to the supplied chromosome.")
  println("  --truncateTo n        limit the number of records retrieved for each reads file.")
  println("  --stageLocal          download the raw reads files to the path specifed in the RNASeq runsheet.  Set to false to disable raw read download and processing.")
  exit 0
  }

println "PARAMS: $params"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
// Get all params sourced data into channels
// Set up channel containing glds accession number
if ( params.gldsAccession ) {ch_glds_accession = Channel.from( params.gldsAccession )} else { exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"}
if ( !params.ensemblVersion ) { exit 1, "Missing Required Parameter: ensemblVersion. Example for setting on CLI: --ensemblVersion 96" }

/**************************************************
* DEBUG WARNING  **********************************
**************************************************/
if ( params.limitSamplesTo || params.truncateTo) {
  println("${c_back_bright_red}WARNING WARNING: DEBUG OPTIONS ENABLED!")
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  println("WARNING WARNING: DEBUG OPTIONS ENABLED!${c_reset}")
} else {
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
}

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/
if ( params.stageLocal && params.truncateTo ) {
  // download truncated raw reads
  println("${c_bright_green}Staging truncated raw reads for ${params.gldsAccession}${c_reset}")
} else if ( params.stageLocal && !params.truncateTo ) {
  // download full raw reads
  println("${c_bright_green}Staging raw reads for ${params.gldsAccession}${c_reset}")
} else {
  // maybe print some nice data from the samplesheet
  println("${c_bright_green}No Staging of raw reads.  Only getting Metadata for ${params.gldsAccession}${c_reset}")
}

include { staging as STAGING } from './stage_analysis.nf'

workflow {
	main:
    STAGING( ch_glds_accession )
    if ( params.stageLocal ) {
      STAGING.out.runsheet
      STAGING.out.raw_reads | set { raw_reads_ch }
      // meta only for dataset specific processes that don't use samples
      // e.g. downloading correct reference genome base on organism
      STAGING.out.raw_reads | take(1) | map{it -> it[0]} | set { meta_ch }
      STAGING.out.isa | set { isa_ch }

      println (meta_ch)

      raw_reads_ch | RAW_FASTQC //| view {"POST_FASTQC: $it"}

      RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
                            //| view {"POST_MAP: $it"}
                            | flatten
                            | unique
                            | collect
                            //| view {"PRE_RAW_MULTIQC: $it"}
                            | RAW_MULTIQC

      raw_reads_ch |  TRIMGALORE

      TRIMGALORE.out.reads | TRIMMED_FASTQC

      TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } \
                                | flatten \
                                | unique \
                                | collect \
                                | TRIMMED_MULTIQC

      meta_ch | DOWNLOAD_GENOME_ANNOTATIONS | set { genome_annotations_pre_subsample }

      // SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
      if ( params.genomeSubsample ) {
        SUBSAMPLE_GENOME( genome_annotations_pre_subsample, meta_ch )
        SUBSAMPLE_GENOME.out.build | first | set { genome_annotations_pre_ercc }
      } else {
        genome_annotations_pre_subsample | first | set { genome_annotations_pre_ercc }
      }

      // ERCC STEP : ADD ERCC Fasta and GTF to genome files
      CONCAT_ERCC( genome_annotations_pre_ercc, DOWNLOAD_ERCC(), meta_ch )
      .ifEmpty { genome_annotations_pre_ercc.value }  | set { genome_annotations }
      meta_ch | view
      BUILD_STAR( genome_annotations, meta_ch)

      TRIMGALORE.out.reads | combine( BUILD_STAR.out.build ) | ALIGN_STAR

      BUILD_RSEM( genome_annotations, meta_ch)

      ALIGN_STAR.out | combine( BUILD_RSEM.out.build ) | set { aligned_ch }
      aligned_ch | COUNT_ALIGNED

      ALIGN_STAR.out | map { it -> it[1] } | collect | ALIGN_MULTIQC

      COUNT_ALIGNED.out | map { it[0].id }
                        | collectFile(name: "samples.txt", newLine: true)
                        | set { samples_ch }

      COUNT_ALIGNED.out | map { it[1] } | collect | set { rsem_ch }

      // TODO: Reintegrate QUANTIFY_GENES( samples_ch, rsem_ch )

      organism_ch = channel.fromPath( params.organismCSV )

      DGE_BY_DESEQ2( isa_ch, organism_ch, rsem_ch, meta_ch  )


      // Software Version Capturing
      ch_software_versions = Channel.empty()
      RAW_FASTQC.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      RAW_MULTIQC.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      TRIMGALORE.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      TRIMMED_FASTQC.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      TRIMMED_MULTIQC.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      BUILD_STAR.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      BUILD_RSEM.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      DGE_BY_DESEQ2.out.version.ifEmpty(null) | mix(ch_software_versions) | set{ch_software_versions}
      ch_software_versions | collectFile(name: "${params.gldsAccession}/software_versions.txt", newLine: true)

      // VV processes
      if ( !params.skipVV ) {
        ch_vv_log_00 =  Channel.fromPath("nextflow_vv_log.tsv")
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
    }
}

workflow.onComplete {
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if ( workflow.success ) {
      println "Raw and Processed data location: ${ params.gldsAccession }"
      println "V&V logs location: ${ params.gldsAccession }/VV_Log${c_reset}"
    }
}
