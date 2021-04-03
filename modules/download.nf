/*
 * Dowload raw paired end reads from Genelab website, assumes a constant format on
 *   Genelab side
 */

process DOWNLOAD_RAW_READS {
  label 'networkBound'
  storeDir "${params.publishDirPath}/00-RawData/Fastq"

  input:
    val(sample)
  output:
    tuple val(sample), path("${sample}_R?_raw.fastq.gz"), emit: raw_reads

  script:
    """
    wget --no-check-certificate --quiet \
    -O ${sample}_R1_raw.fastq.gz \
    ${params.GLDS_URL_PREFIX}${sample}_R1_raw.fastq.gz${params.GLDS_URL_SUFFIX}
    """
    if ( params.pairedEnd ) {
      """
      wget --no-check-certificate --quiet \
      -O ${sample}_R2_raw.fastq.gz \
      ${params.GLDS_URL_PREFIX}${sample}_R2_raw.fastq.gz${params.GLDS_URL_SUFFIX}
      """
    }
}

/*
 * Download and decompress genome and annotation files
 */

process DOWNLOAD_GENOME_ANNOTATIONS {\
  conda "${baseDir}/envs/download_tools.yml"
  tag "Organism: ${ meta.organism_sci }  Ensembl Version: ${params.ensemblVersion}"
  label 'networkBound'
  storeDir "${params.storeDirPath}/ensembl/${params.ensemblVersion}/${ meta.organism_sci }"

  input:
    val(meta)
  output:
    tuple path("*.dna.toplevel.fa"), path("*.${ params.ensemblVersion }.gtf")
  script:
    """
    retrieve_references.py --ensembl_version ${ params.ensemblVersion } \
                           --organism        ${ meta.organism_sci}

    # decompress files
    gunzip *.fa.gz
    gunzip *.gtf.gz
    """
}

process DOWNLOAD_ERCC {
  label 'networkBound'
  storeDir "${params.storeDirPath}/thermofisher/ERCC"

  input:
  output:
    tuple path("ERCC92.fa"), path("ERCC92.gtf")
  script:
    """
    wget --no-check-certificate --quiet \
    -O ERCC92.zip \
    https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip  \
    && \
    unzip ERCC92.zip
    """
}
// DEPRECATED
process DOWNLOAD_ISA {
  conda "${baseDir}/envs/download_tools.yml"
  publishDir "${params.publishDirPath}/${ params.metaDataPath }"

  input:
    val(accession)

  output:
    path("*.zip")

  script:
    """
    retrieve_isa_from_genelab.py --accession ${accession} --alternate_url
    """

}
