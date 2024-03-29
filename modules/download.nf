/*
 * Processes for downloading reference files
 */

process DOWNLOAD_GENOME_ANNOTATIONS {
  // Download and decompress genome and annotation files
  tag "Organism: ${ organism_sci }  Ensembl Version: ${params.ensemblVersion}"
  label 'networkBound'
  storeDir "${params.referenceStorePath}/ensembl/${params.ensemblVersion}/${ organism_sci }"

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
    gunzip *.fa.gz || true # this is to accomodate dummy marker file, D.N.E., which results in no decompression needed.
    gunzip *.gtf.gz
    """
}

process DOWNLOAD_ERCC {
  label 'networkBound'
  storeDir "${params.referenceStorePath}/ERCC_thermofisher"

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
