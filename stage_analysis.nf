/*
* Workflow that accepts a GLDS accession and generates the following:
* 1. Download ISA.zip and generates RNASeq Samplesheet
* 2a. Downloads Raw Reads
* 2b. Downloads Truncated Raw Reads (Useful for testing with limited resources)
*/
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";
// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2


// Import process from separate module file
include { RNASEQ_SAMPLESHEET_FROM_GLDS as GENERATE_SAMPLESHEET;
          STAGE_RAW_READS;
          get_samplesheet_paths } from'./modules/genelab.nf'

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ RNASeq Concensus Pipeline: Staging Workflow: $workflow.manifest.version ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("usage: nextflow run J-81/Nextflow_RCP -main-script stage_analysis.nf [-h] [--limitSamplesTo n] [--truncateTo n] [--stageLocal] --gldsAccession GLDS-000")
  println()
  println("required arguments:")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession number to stage raw reads for the RNASeq Concensus Pipeline")
  println("optional arguments:")
  println("  --help                show this help message and exit")
  println("  --limitSamplesTo n    limit the number of samples staged to a number.")
  println("  --truncateTo n        limit the number of records retrieved for each reads file.")
  println("  --stageLocal          download the raw reads files to the path specifed in the RNASeq samplesheet.")
  exit 0
  }
/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
// Get all params sourced data into channels
// Set up channel containing glds accession number
if (params.gldsAccession) {ch_glds_accession = Channel.from( params.gldsAccession )} else { exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"}

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
  println("${c_bright_green}No Staging.  Only getting Metadata for ${params.gldsAccession}${c_reset}")
}

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow {
  main:
    sample_limit = params.limitSamplesTo ? params.limitSamplesTo : -1 // -1 in take means no limit

    ch_glds_accession | GENERATE_SAMPLESHEET

    GENERATE_SAMPLESHEET.out.samplesheet | splitCsv(header: true)
                                | map{ row -> get_samplesheet_paths(row) }
                                | take( sample_limit )
                                | set{ ch_samples }

    if ( params.stageLocal && params.truncateTo ) {
      // download truncated raw reads
      // download full raw reads
      ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]]}
                 | branch {
                   paired: it.size() == 3
                   single: it.size() == 2
                 }
                 | set{ ch_raw_read_pointers}

       // PAIRED END
       // Only difference is the splitFastq arg 'pe'
      ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: params.truncateTo, by: params.truncateTo, file: true)
                                  | map { it -> [ it[0], [ it[1], it[2] ] ]}
                                  //| view { it -> "TRUNCATED PAIRED READS ($params.truncateTo): $it[0]"}
                                  | set { ch_raw_reads }
       // SINGLE END
       // Only difference is the splitFastq arg 'pe'
      ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: params.truncateTo, by: params.truncateTo, file: true)
                                  | map { it -> [ it[0], [ it[1] ] ]}
                                  //| view { it -> "TRUNCATED SINGLE READS ($params.truncateTo): $it[0]"}
                                  | mix( ch_raw_reads )
                                  | set { ch_raw_reads }

      // Moves the truncated files to expected raw read locations as per samplesheet
      ch_raw_reads | STAGE_RAW_READS

    } else if ( params.stageLocal && !params.truncateTo ) {
      // download full raw reads
      ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                 | set { ch_raw_reads }

      // Download the raw reads and publish them to expected raw read locations as per samplesheet
    ch_raw_reads | STAGE_RAW_READS

    } else {
      // download nothing, end of workflow
      // maybe print some nice data from the samplesheet

    }
}
