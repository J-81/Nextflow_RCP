/*
 * Dowload raw paired end reads from Genelab website, assumes a constant format on
 *   Genelab side
 */

process DOWNLOAD_RAW_READS {
  label 'networkBound'
  publishDir "${params.publishDirPath}/00-RawData/Fastq"

  input:
    val(sample)
  output:
    tuple val(sample), path("${sample}_R?_raw.fastq.gz"), emit: raw_reads

  stub:
    if ( params.pairedEnd ) {
      """
      touch ${sample}_R1_raw.fastq.gz
      touch ${sample}_R2_raw.fastq.gz
      """
    } else {
      """
      touch ${sample}_R1_raw.fastq.gz
      """
    }

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
  storeDir "${params.storeDirPath}/ensembl"

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
