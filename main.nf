nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { FASTQC as RAW_FASTQC } from './modules/quality.nf' addParams(PublishTo: "00-RawData/FastQC_Reports")
include { FASTQC as TRIMMED_FASTQC } from './modules/quality.nf' addParams(PublishTo: "01-TG_Preproc/FastQC_Reports")
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "00-RawData/FastQC_Reports", MQCLabel:"raw")
include { MULTIQC as TRIMMED_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "01-TG_Preproc/FastQC_Reports", MQCLabel:"trimmed")
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "01-TG_Preproc/Trimming_Reports", MQCLabel:"trimming")
include { MULTIQC as ALIGN_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "02-STAR_Alignment", MQCLabel:"align")
include { MULTIQC as COUNT_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "03-RSEM_Counts", MQCLabel:"count")
include { MULTIQC as ALL_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "MULTIQC_ALL", MQCLabel:"all")
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_STAR;
          ALIGN_STAR;
          BUILD_RSEM;
          COUNT_ALIGNED;
          SUBSAMPLE_GENOME;
          CONCAT_ERCC;
          QUANTIFY_GENES } from './modules/genome.nf'
include { DGE_BY_DESEQ2 } from './modules/dge.nf'
include { VV_RAW_READS;
          VV_TRIMMED_READS;
          VV_RAW_READS_MULTIQC;
          VV_TRIMMED_READS_MULTIQC;
          VV_STAR_ALIGNMENTS;
          VV_RSEQC;
          VV_RSEM_COUNTS;
          VV_DESEQ2_ANALYSIS } from './modules/vv.nf' addParams( RootDirForVV: "${workflow.launchDir}/${ params.outputDir }")
include { GET_MAX_READ_LENGTH } from './modules/fastqc.nf'
include { POST_PROCESSING } from './modules/genelab.nf'

/**************************************************
* HELP MENU  **************************************
**************************************************/
allowed_ref_order = ['toplevel','primary_assemblyELSEtoplevel']
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ RNASeq Concensus Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("usage: nextflow run J-81/Nextflow_RCP -r $workflow.revision --gldsAccession GLDS-000 --ensemblVersion 96  [--skipVV] [--outputDir] [--stageLocal] [--limitSamplesTo n] [--truncateTo n] [--genomeSubsample n]")
  println()
  println("required arguments:")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession id to process through the RNASeq Concensus Pipeline.")
  println("  --ensemblVersion n    the ensembl Version to use for the reference genome.")
  println("optional arguments:")
  println("  --help                show this help message and exit")
  println("  --skipVV              skip automated V&V processes. Default: false")
  println("  --outputDir           directory to save staged raw files and processed files. Default: <launch directory>")
  println("  --limitSamplesTo n    limit the number of samples staged to a number.")
  println("  --genomeSubsample n   subsamples genome fasta and gtf files to the supplied chromosome.")
  println("  --truncateTo n        limit number of reads downloaded and processed to *n* reads , for paired end limits number of reverse and forward read files to *n* reads each.")
  println("  --force_single_end    forces analysis to use single end processing.  For paired end datasets, this means only R1 is used.  For single end studies, this should have no effect.")
  println("  --stageLocal          download the raw reads files for the supplied GLDS accession id.  Set to false to disable raw read download and processing.  Default: true")
  println("  --ref_order           specifies the reference to use from ensembl.  Allowed values:  ['toplevel','primary_assemblyELSEtoplevel']. 'toplevel' : use toplevel.  'primary_assemblyELSEtoplevel' : use primary assembly, but use toplevel if primary assembly doesn't exist. Default: 'primary_assemblyELSEtoplevel'")  
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

if ( !allowed_ref_order.contains(params.ref_order ) ) { exit 1, "Invalid ref_order param.  Must be either 'toplevel' or 'primary_assembly,toplevel'" }
if ( !params.outputDir ) {  params.outputDir = "$workflow.launchDir" }

ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")

/**************************************************
* DEBUG WARNING  **********************************
**************************************************/
if ( params.limitSamplesTo || params.truncateTo || params.force_single_end) {
  println("${c_back_bright_red}WARNING WARNING: DEBUG OPTIONS ENABLED!")
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  params.force_single_end ? println("Forcing analysis to used only forward reads if paired end (i.e. as though single ended") : println("No forcing single end analysis")
  println("WARNING WARNING: DEBUG OPTIONS ENABLED!${c_reset}")
} else {
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  params.force_single_end ? println("Forcing analysis to used only forward reads if paired end (i.e. as though single ended") : println("No forcing single end analysis")
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
include { references as REFERENCES } from './references.nf'
include { strandedness as STRANDEDNESS } from './strandedness.nf'

workflow {
	main:
    STAGING( ch_glds_accession )
    if ( params.stageLocal ) {
      STAGING.out.runsheet
      STAGING.out.raw_reads | set { raw_reads_ch }
      // meta only for dataset specific processes that don't use samples
      // e.g. downloading correct reference genome base on organism
      STAGING.out.raw_reads | take(1) | map{it -> it[0]} | set { meta_ch }
      STAGING.out.raw_reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
                            | set { samples_ch }
      STAGING.out.isa | set { isa_ch }

      meta_ch | view

      raw_reads_ch | RAW_FASTQC //| view {"POST_FASTQC: $it"}

      RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
                            //| view {"POST_MAP: $it"}
                            | flatten
                            | unique
                            | collect
                            //| view {"PRE_RAW_MULTIQC: $it"}
                            | set { raw_mqc_ch }

      RAW_FASTQC.out.fastqc | map { it -> [ it[2] ] }
                            | flatten
                            | GET_MAX_READ_LENGTH

      GET_MAX_READ_LENGTH.out.length | max { it.toInteger() }
                                     | set { max_read_length_ch }

      raw_reads_ch |  TRIMGALORE

      TRIMGALORE.out.reads | TRIMMED_FASTQC

      TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } \
                                | flatten \
                                | unique \
                                | collect \
                                | set { trim_mqc_ch }

      REFERENCES( meta_ch | map { it.organism_sci }, meta_ch | map { it.has_ercc } )
      REFERENCES.out.genome_annotations | set { genome_annotations }      

      BUILD_STAR( genome_annotations, meta_ch, max_read_length_ch)

      TRIMGALORE.out.reads | combine( BUILD_STAR.out.build ) | ALIGN_STAR

      STRANDEDNESS ( ALIGN_STAR.out.bam_by_coord, REFERENCES.out.genome_bed, samples_ch ) 
      STRANDEDNESS.out.strandedness | map { it.text.split(":")[0] } | set { strandedness_ch }

      BUILD_RSEM( genome_annotations, meta_ch)

      ALIGN_STAR.out.alignments | combine( BUILD_RSEM.out.build ) | set { aligned_ch }
      COUNT_ALIGNED( aligned_ch, strandedness_ch )

      ALIGN_STAR.out.alignments | map { it -> it[1] } 
                                | collect 
                                | set { align_mqc_ch }


      COUNT_ALIGNED.out.counts | map { it[1] } | collect | set { rsem_ch }

      QUANTIFY_GENES( samples_ch, rsem_ch )

      organism_ch = channel.fromPath( params.organismCSV )

      DGE_BY_DESEQ2( STAGING.out.runsheet, organism_ch, rsem_ch, meta_ch  )


      // ALL MULTIQC
      RAW_MULTIQC( samples_ch, raw_mqc_ch, ch_multiqc_config  )
      TRIMMED_MULTIQC( samples_ch, trim_mqc_ch, ch_multiqc_config ) // refering to the trimmed reads
      TRIM_MULTIQC( samples_ch, TRIMGALORE.out.reports, ch_multiqc_config ) // refering to the trimming process
      ALIGN_MULTIQC( samples_ch, align_mqc_ch, ch_multiqc_config )
      COUNT_MULTIQC( samples_ch, rsem_ch, ch_multiqc_config )
      raw_mqc_ch | concat( trim_mqc_ch ) 
                 | concat( ALIGN_STAR.out.alignment_logs ) 
                 | concat( STRANDEDNESS.out.rseqc_logs )
                 | concat( rsem_ch )
                 | concat( TRIMGALORE.out.reports )
                 | collect | set { all_mqc_ch }
      ALL_MULTIQC( samples_ch, all_mqc_ch, ch_multiqc_config )

      // Software Version Capturing
      nf_version = "Nextflow Version:".concat("${nextflow.version}\n<><><>\n")
      ch_nextflow_version = Channel.value(nf_version)
      ch_software_versions = Channel.empty()
      RAW_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
      RAW_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
      TRIMGALORE.out.version | mix(ch_software_versions) | set{ch_software_versions}
      TRIMMED_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
      TRIMMED_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
      ALIGN_STAR.out.version | mix(ch_software_versions) | set{ch_software_versions}
      COUNT_ALIGNED.out.version | mix(ch_software_versions) | set{ch_software_versions}
      DGE_BY_DESEQ2.out.version | mix(ch_software_versions) | set{ch_software_versions}
      STRANDEDNESS.out.versions | mix(ch_software_versions) | set{ch_software_versions}
      ch_software_versions | map { it.text + "\n<><><>\n"}
                           | unique
                           | mix(ch_nextflow_version)
                           | collectFile(name: "software_versions.txt", newLine: true, cache: false)
			   | set{ch_final_software_versions}


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

        VV_STAR_ALIGNMENTS( ALIGN_STAR.out.alignments | map{ it -> it[1..it.size()-1] } | collect, // map use here: removes val(meta) from tuple
                            ch_vv_log_04 ) | set { ch_vv_log_05 }
        
        //VV_RSEQC( STRANDEDNESS.out.infer_expt_mqc, 
        //                    ch_vv_log_05 ) | set { ch_vv_log_06 }

        ch_vv_log_05 | set { ch_vv_log_06 } // temporary skip

        VV_RSEM_COUNTS( COUNT_ALIGNED.out.counts | map{ it -> it[1..it.size()-1] } | flatten | collect, // map use here: removes val(meta) from tuple
                        ch_vv_log_06 ) | set { ch_vv_log_07 }

        VV_DESEQ2_ANALYSIS( DGE_BY_DESEQ2.out.dge | map{ it -> it[1..it.size()-1] } | collect, // map use here: removes val(meta) from tuple
                            ch_vv_log_07 ) | set { ch_vv_log_08 }
        
        // GeneLab post processing
        POST_PROCESSING(STAGING.out.runsheet, ch_final_software_versions, ch_vv_log_08, STAGING.out.metasheet) // Penultimate process when V&V enabled is the last V&V process
      } else {
        POST_PROCESSING(STAGING.out.runsheet, ch_final_software_versions, Channel.value("NO VV, last output is software versions"), STAGING.out.metasheet)
      }
    }
}

workflow.onComplete {
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if ( workflow.success ) {
      println "Raw and Processed data location: ${ params.outputDir }/${ params.gldsAccession }"
      println "V&V logs location: ${ params.outputDir }/${ params.gldsAccession }/VV_Log"
      println "Pipeline tracing/visualization files location: ${ params.tracedir }/${ params.gldsAccession }${c_reset}"
    }
}
