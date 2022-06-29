/* Processes dealing with retrieving data from GeneLab
*/
// TODO Migrate these CLI args to updated dp tools package API
process RNASEQ_RUNSHEET_FROM_GLDS {
  // Downloads isazip and creates run sheets using GeneLab API
  tag "${ glds_accession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(glds_accession)

  output:
    path("${ glds_accession }_bulkRNASeq_v*.csv"), emit: runsheet
    path("*.zip"), emit: isazip

  script:
    """
    dpt-get-isa-archive --accession ${ glds_accession }\
      --alternate-url

    dpt-isa-to-runsheet --accession ${ glds_accession } \
      --config-type bulkRNASeq --isa-archive *.zip
    """
}


process STAGE_RAW_READS {
  // Stages the raw reads into appropriate publish directory
  tag "${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.raw_reads_root_dir }/Fastq",
    mode: params.publish_dir_mode

  input:
    tuple val(meta), path("?.gz")

  output:
    tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

  script:
    if ( meta.paired_end ) {
      """
      cp -L 1.gz ${meta.id}_R1_raw.fastq.gz
      cp -L 2.gz ${meta.id}_R2_raw.fastq.gz
      """
    } else {
      """
      cp -L 1.gz  ${meta.id}_raw.fastq.gz
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
    path(runsheet)

  output:
    path("${ params.gldsAccession }_metadata_table.txt"), emit: metasheet

  script:
    """
    create_table_v2.py --accession ${ params.gldsAccession }  \
                       --isa-zip isa.zip \
                       --output-dir . \
                       --runsheet ${ runsheet }
    """
}

process POST_PROCESSING {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path("runsheet.csv")
    val(LAST_PROCESS_MARKER) // Unused in task, but used in workflow definition to ensure this process is last regardless of whether V&V is used
    val(LAST_PROCESS_MARKER_2) // Unused in task, but used in workflow definition to ensure this process follows metasheet generation

  output:
    path("updated_curation_tables") // directory containing extended ISA tables
    path("*md5sum*")

  script:
    root_out_dir = "${ workflow.launchDir}/${ params.outputDir }/${ params.gldsAccession }"
    """
    generate_md5sum_files.py --root-path ${ root_out_dir } --accession ${ params.gldsAccession }
    update_curation_table.py --root-path ${ root_out_dir } --accession ${ params.gldsAccession }
    """
}

process SOFTWARE_VERSIONS {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path("software_versions.txt")

  output:
    path("software_versions.md")

  script:
    """
    format_software_versions.py software_versions.txt
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

    def PRIMARY_KEYS = ["mus_musculus":"ENSEMBL",
                        "danio_rerio":"ENSEMBL",
                        "rattus_norvegicus":"ENSEMBL",
                        "homo_sapiens":"ENSEMBL",
                        "drosophila_melanogaster":"ENSEMBL",
                        "caenorhabditis_elegans":"ENSEMBL",
                        "arabidopsis_thaliana":"TAIR"]

    def meta = [:]
    meta.id                         = row["Sample Name"]
    meta.organism_sci               = row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.primary_keytype            = PRIMARY_KEYS[meta.organism_sci]
    meta.paired_end                 = row.paired_end.toBoolean()
    meta.has_ercc                   = row.has_ERCC.toBoolean()

    def array = []
    def raw_reads = []
    raw_reads.add(file(row.read1_path))
    if (meta.paired_end) {
        raw_reads.add(file(row.read2_path))
      }
    array = [meta, raw_reads]
    return array
}
