/*
 * Processes for downloading reference files
 */

process DOWNLOAD_GENOME_ANNOTATIONS {
  // Download and decompress genome and annotation files
  tag "storeDir: ${ task.storeDir } Target(s): *.dna.${ params.ref_target }.fa, *.${ params.ensemblVersion }.gtf Note: exact target name undefined at task start, hence '*' here"
  label 'networkBound'
  storeDir "${params.referenceStorePath}/ensembl/${params.ensemblVersion}/${ organism_sci }"

  input:
    val(organism_sci)

  output:
    tuple path("*.dna.${ params.ref_target }.fa"), path("*.${ params.ensemblVersion }.gtf")

  script:
    """
    # ensure proper permissions for generated file
    umask 0022

    retrieve_references.py --ensembl_version ${ params.ensemblVersion } \
                           --organism        ${ organism_sci} \
                           --target          ${ params.ref_target }



    # decompress files
    gunzip *.fa.gz || true # this is to accomodate dummy marker file, D.N.E., which results in no decompression needed.
    gunzip *.gtf.gz
    """
}

process DOWNLOAD_ERCC {
  tag "storeDir: ${ task.storeDir } Target(s): ERCC92.fa, ERCC92.gtf"
  label 'networkBound'
  storeDir "${params.referenceStorePath}/ERCC_thermofisher"

  input:

  output:
    tuple path("ERCC92.fa"), path("ERCC92.gtf")

  script:
    """
    # ensure proper permissions for generated file
    umask 0022

    wget --no-check-certificate --quiet \
    -O ERCC92.zip \
    https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip  \
    && \
    unzip ERCC92.zip
    """
}
