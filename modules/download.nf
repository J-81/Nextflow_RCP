/*
 * Processes for downloading reference files
 */

process DOWNLOAD_GENOME_ANNOTATIONS {
  // Download and decompress genome and annotation files
  tag "Organism: ${ organism_sci }  Ensembl Version: ${params.ensemblVersion}"
  label 'networkBound'
  storeDir "${params.storeDirPath}/ensembl/${params.ensemblVersion}/${ organism_sci }"
  errorStrategy "${ params._has_fallback }" ? 'ignore' : 'terminate'

  input:
    val(organism_sci)

  output:
    tuple path("*.dna.${ params.ref_target }.fa"), path("*.${ params.ensemblVersion }.gtf")

  script:
    """
    retrieve_references.py --ensembl_version ${ params.ensemblVersion } \
                           --organism        ${ organism_sci} \
                           --target          ${ params.ref_target }



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
