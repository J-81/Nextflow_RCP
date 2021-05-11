/*
* Workflow that accepts a GLDS accession and generates the following:
* 1. Download ISA.zip and generates RNASeq Samplesheet
* 2a. Downloads Raw Reads
* 2b. Downloads Truncated Raw Reads (Useful for testing with limited resources)
*/

// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2


// Import process from separate module file
include { RNASEQ_RUNSHEET_FROM_GLDS as GENERATE_RUNSHEET;
          GENERATE_METASHEET;
          STAGE_RAW_READS;
          get_runsheet_paths } from'./modules/genelab.nf'

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow staging{
  take:
    ch_glds_accession
  main:
    sample_limit = params.limitSamplesTo ? params.limitSamplesTo : -1 // -1 in take means no limit

    ch_glds_accession | GENERATE_RUNSHEET

    GENERATE_RUNSHEET.out.runsheet | splitCsv(header: true)
                                   | map{ row -> get_runsheet_paths(row) }
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
      // Don't download any raw reads
    }

    GENERATE_RUNSHEET.out.isazip | GENERATE_METASHEET

    emit:
      raw_reads = params.stageLocal ? STAGE_RAW_READS.out : null
      isa = GENERATE_RUNSHEET.out.isazip
      runsheet = GENERATE_RUNSHEET.out.runsheet
}
