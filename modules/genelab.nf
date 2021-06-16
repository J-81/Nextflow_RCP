/* Processes dealing with retrieving data from GeneLab
*/
process RNASEQ_RUNSHEET_FROM_GLDS {
  // Downloads isazip and creates run sheets using GeneLab API
  tag "${ glds_accession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(glds_accession)

  output:
    path("AST_autogen_*_${ glds_accession }_RNASeq_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isazip

  script:
    """
    retrieve_isa_from_genelab.py --accession ${ glds_accession }\
                                 --to-RNASeq-runsheet
    """
}


process STAGE_RAW_READS {
  // Stages the raw reads into appropriate publish directory
  tag "${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ meta.raw_read_root_dir }",
    mode: params.publish_dir_mode

  input:
    tuple val(meta), path("?.gz")

  output:
    tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

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


process GENERATE_METASHEET {
  // Generates a metadata table, not used in further processing
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    path("isa.zip")

  output:
    path("${ params.gldsAccession }_metadata_table.txt")

  script:
    """
    create_table_v2.py --accession ${ params.gldsAccession }  \
                       --isa-zip isa.zip \
                       --output-dir .
    """
}

// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Original Function Credit: Dr. Harshil Patel
// Function to get list of [ meta, [ fastq_1_path, fastq_2_path ] ]
def get_runsheet_paths(LinkedHashMap row) {
    def ORGANISMS = ["mus_musculus":"MOUSE",
                     "danio_rerio":"ZEBRAFISH",
                     "rattus_norvegicus":"RAT",
                     "homo_sapiens":"HUMAN",
                     "drosophila_melanogaster":"FLY",
                     "caenorhabditis_elegans":"WORM",
                     "arabidopsis_thaliana":"ARABIDOPSIS"]

    def meta = [:]
    meta.id                         = row.sample_name
    meta.organism_sci               = row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.read_length_R1             = row.read_length_R1.toInteger()
    meta.paired_end                 = row.paired_end.toBoolean()
    meta.has_ercc                   = row.has_ERCC.toBoolean()
    meta.raw_read1                  = new File(row.raw_read1) //points to file
    meta.raw_read_root_dir          = new File(row.raw_read1).parent //points to file, parent is directory to store
    meta.raw_read_fastQC            = new File(row.raw_read_fastQC) //points to directory
    meta.raw_read_multiqc           = new File(row.raw_read_multiqc).parent //points to file, parent is directory to store Raw Read MultiQC output to
    meta.trimmed_read_root_dir      = new File(row.trimmed_read1).parent //points to file, parent is directory to store
    meta.trimmed_read1              = new File(row.trimmed_read1) //points to file
    meta.trimmed_read_fastQC        = new File(row.trimmed_read_fastQC) //points to directory
    meta.trimmed_read_multiqc       = new File(row.trimmed_read_multiqc).parent //points to file, parent is directory to store Trimmed Read MultiQC output to
    meta.STAR_Alignment_dir         = new File(row.STAR_Alignment) //points to sample directory
    meta.STAR_Alignment_root_dir    = new File(row.STAR_Alignment).parent //points to sample directory
    meta.RSEM_Counts_dir            = new File(row.RSEM_Counts) //points to sample directory
    meta.RSEM_Counts_root_dir       = new File(row.RSEM_Counts).parent //points to sample directory
    meta.DESeq2_NormCount           = new File(row.DESeq2_NormCount) //points to dataset directory, same for all samples
    meta.DESeq2_DGE                 = new File(row.DESeq2_DGE) //points to dataset directory, same for all samples

    def array = []
    def raw_reads = []
    meta.stage1 = file("${ params.gldsAccession }") / file(row.raw_read1).name
    raw_reads.add(file(row.read1_url))
    if (meta.paired_end) {
        meta.stage2                         = file("${ params.gldsAccession }") / file(row.raw_read2).name
        meta.raw_read2                      = new File(row.raw_read2) //points to file
        meta.trimmed_read2                  = new File(row.trimmed_read2) //points to file
        meta.read_length_R2                 = row.read_length_R2.toInteger()
        raw_reads.add(file(row.read2_url))
      }
    array = [meta, raw_reads]
    return array
}
