/* Processes dealing with retrieving data from GeneLab
*/

process RNASEQ_SAMPLESHEET_FROM_GLDS {
  conda "${baseDir}/envs/AST"
  tag "${ glds_accession }"
  storeDir "${params.gldsAccession}/Metadata"

  input:
    val(glds_accession)

  output:
    path("AST_autogen_${ glds_accession }_RNASeq_samplesheet.csv"), emit: samplesheet
    path("*.zip"), emit: isazip

  script:
    """
    retrieve_isa_from_genelab.py --accession ${ glds_accession }\
                                 --alternate_url\
                                 --to_RNASeq_samplesheet
    """

}


process STAGE_RAW_READS {
  tag "${ meta.id }"
  storeDir "${ params.gldsAccession }/00-RawData/Fastq"

  input:
    tuple val(meta), path("?.gz")

  output:
    tuple val(meta), path("${meta.id}*.gz")

  script:
    if ( meta.paired_end ) {
      """
      cp -L 1.gz ${meta.stage1.name}
      cp -L 2.gz ${meta.stage2.name}
      """
    } else {
      """
      cp -L 1.gz ${meta.stage1.name}
      """
    }
}

// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Original Function Credit: Dr. Harshil Patel
// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample_name
    meta.paired_end   = row.paired_end.toBoolean()
    meta.has_ercc     = row.has_ERCC.toBoolean()
    //meta.strandedness = row.strandedness

    def array = []
    def raw_reads = []
    meta.stage1 = file("${ params.gldsAccession }") / file(row.raw_read1).name
    raw_reads.add(file(row.read1_url))
    if (meta.paired_end) {
        meta.stage2 = file("${ params.gldsAccession }") / file(row.raw_read2).name
        raw_reads.add(file(row.read2_url))
      }
    array = [meta, raw_reads]
    return array
}
