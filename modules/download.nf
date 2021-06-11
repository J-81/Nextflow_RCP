/*
 * Processes for downloading reference files
 */

process DOWNLOAD_GENOME_ANNOTATIONS {
  // Download and decompress genome and annotation files
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
