/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  conda "${baseDir}/envs/fastqc.yml"
  cpus { read.size() } // number of read files to process
  //storeDir "${params.storeDirPath}/fastqc"

  input:
    tuple val(sample), path(read)
  output:
    tuple val(sample), path("*_fastqc.html"), path("*_fastqc.zip")

  stub:
    if ( params.pairedEnd ) {
      """
      # Forward Reads
      touch ${sample}_R1_raw_fastqc.html
      touch ${sample}_R1_raw_fastqc.zip

      # Reverse Reads
      touch ${sample}_R2_raw_fastqc.html
      touch ${sample}_R2_raw_fastqc.zip
      """
    } else {
      """
      # Single Reads (labeled like forward reads)
      touch ${sample}_R1_raw_fastqc.html
      touch ${sample}_R1_raw_fastqc.zip
      """
    }

  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $read
    """

}

process MULTIQC {
  label "fastLocal"
  conda "${baseDir}/envs/multiqc.yml"
  publishDir "${params.publishDirPath}/multiQC/${params.multiQCLabel}"

  input:
    path(fastqc) // any number of fastqc files
  output:
    path("${params.multiQCLabel}_multiqc_report.html")

  stub:
    """
    touch "${params.multiQCLabel}_multiqc_report.html"
    """

  script:
    """
    multiqc -n ${params.multiQCLabel}_multiqc_report .
    """

}

process TRIMGALORE {
  conda "${baseDir}/envs/trim_galore.yml"
  cpus 4

  input:
    tuple val(sample), path(forward_read), path(reverse_read)
  output:
    tuple val(sample), path("${ forward_read.simpleName }_val_1.fq.gz"), path("${ reverse_read.simpleName }_val_2.fq.gz"), emit: reads
    tuple val(sample), path("${ forward_read }_trimming_report.txt"), path("${ reverse_read }_trimming_report.txt"), emit: trim_reports

  stub:
  // TODO: add support for single end studies
    """
    touch "${ forward_read.simpleName }_val_1.fq.gz"
    touch "${ reverse_read.simpleName }_val_2.fq.gz"

    touch "${ forward_read }_trimming_report.txt"
    touch "${ reverse_read }_trimming_report.txt"
    """

  script:
    /*
     * comments -> --ilumina # if adapters are not illumina, replace with adapters
     *   --paired  # only for PE studies, # if SE use only single read file
     */
    """
    trim_galore --gzip \
    --cores $task.cpus \
    --illumina \
    --phred33 \
    --paired $forward_read $reverse_read \
    --output_dir .
    """
}
