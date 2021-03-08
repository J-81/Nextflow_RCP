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

process DOWNLOAD_GENOME_ANNOTATIONS {
  label 'networkBound'
  storeDir "${params.storeDirPath}/ensembl/${params.ensembl_version}/${params.organism}"

  input:
  output:
    tuple path("Mus_musculus.GRCm38.dna.toplevel.fa"), path("Mus_musculus.GRCm38.${ params.ensembl_version }.gtf")
  script:
    """
    wget --no-check-certificate --quiet \
    -O Mus_musculus.GRCm38.dna.toplevel.fa.gz \
    ftp://ftp.ensembl.org/pub/release-${ params.ensembl_version }/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz \
    && \
    gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz \
    && \
    wget --no-check-certificate --quiet \
    -O Mus_musculus.GRCm38.${ params.ensembl_version }.gtf.gz \
    ftp://ftp.ensembl.org/pub/release-${ params.ensembl_version }/gtf/mus_musculus/Mus_musculus.GRCm38.${ params.ensembl_version }.gtf.gz \
    && \
    gunzip Mus_musculus.GRCm38.${ params.ensembl_version }.gtf.gz \

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


/* Stub process
TODO: replace with ISA download via api
*/
process DOWNLOAD_ISA {
  publishDir "${params.publishDirPath}/${ params.metaDataPath }"

  input:
  output:
    path("GLDS-${ params.GLDS }_metadata_GLDS-${ params.GLDS }-ISA.zip")

  script:
    """
    cp ${ params.ISAZip } GLDS-${ params.GLDS }_metadata_GLDS-${ params.GLDS }-ISA.zip
    """

}
