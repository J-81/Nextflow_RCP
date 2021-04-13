/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process RAW_FASTQC {
  conda "${baseDir}/envs/fastqc.yml"
  tag "Sample: ${ meta.id }"
  // cpus { read.size() } // BUGGED FOR SINGLE READS: number of read files to process
  publishDir "${ params.gldsAccession }/${ meta.raw_read_fastQC }"

  input:
    tuple val(meta), path(reads)
  output:
    tuple val(meta), path("${ meta.id }*.html"), path("${ meta.id }*.zip"), emit: fastqc
    path("versions.txt"), emit: version

  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $reads

    fastqc -v > versions.txt
    """
}

process TRIMMED_FASTQC {
  conda "${baseDir}/envs/fastqc.yml"
  tag "Sample: ${ meta.id }"
  // cpus { read.size() } // BUGGED FOR SINGLE READS: number of read files to process
  publishDir "${ params.gldsAccession }/${ meta.trimmed_read_fastQC }"

  input:
  tuple val(meta), path(reads)
  output:
  tuple val(meta), path("${ meta.id }*.html"), path("${ meta.id }*.zip"), emit: fastqc
  path("versions.txt"), emit: version

  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $reads

    fastqc -v > versions.txt
    """
}

process RAW_MULTIQC {
  label "fastLocal"
  tag "Dataset: ${ params.gldsAccession }"
  conda "${baseDir}/envs/multiqc.yml"
  publishDir "${ params.gldsAccession }/00-RawData/FastQC_Reports"

  input:
    path("fastqc/*") // any number of fastqc files
  output:
    path("raw_multiqc_report/raw_multiqc.html"), emit: html
    path("raw_multiqc_report/raw_multiqc_data"), emit: data
    path("versions.txt"), emit: version

  script:
    """
    multiqc -o raw_multiqc_report -n raw_multiqc fastqc

    multiqc --version > versions.txt
    """
}

process TRIMMED_MULTIQC {
  label "fastLocal"
  tag "Dataset: ${ params.gldsAccession }"
  conda "${baseDir}/envs/multiqc.yml"
  publishDir "${ params.gldsAccession }/01-TG_Preproc/FastQC_Reports"


  input:
    path("fastqc/*") // any number of fastqc files
  output:
    path("trimmed_multiqc_report/trimmed_multiqc.html"), emit: html
    path("trimmed_multiqc_report/trimmed_multiqc_data"), emit: data
    path("versions.txt"), emit: version

  script:
    """
    multiqc -o trimmed_multiqc_report -n trimmed_multiqc fastqc

    multiqc --version > versions.txt
    """
}

process ALIGN_MULTIQC {
  label "fastLocal"
  //tag "Dataset: ${ params.gldsAccession }"
  conda "${baseDir}/envs/multiqc.yml"
  publishDir "${ params.gldsAccession }/02-STAR_Alignment"

  input:
    path("alignments/*")

  output:
    path("align_multiqc_report/align_multiqc.html"), emit: html
    path("align_multiqc_report/align_multiqc_data"), emit: data

  script:
    """
    multiqc -o align_multiqc_report -n align_multiqc alignments
    """
}

process TRIMGALORE {
  conda "${baseDir}/envs/trim_galore.yml"
  tag "Sample: ${ meta.id }"
  // changing these to storeDir causes an error as pattern is not available for storeDir
  publishDir "${ params.gldsAccession }/01-TG_Preproc/Fastq", pattern: "*trimmed.fastq.gz"
  publishDir "${ params.gldsAccession }/01-TG_Preproc/Trimming_Reports", pattern: "*.txt"

  input:
    tuple val(meta), path(reads)
  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${ meta.id }*.txt"), emit: trim_reports
    path("versions.txt"), emit: version

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
    ${ meta.paired_end ? '--paired' : '' } \
    $reads \
    --output_dir .

    # rename with _trimmed suffix
    ${ meta.paired_end ? \
      "cp ${ meta.id }_R1_raw_val_1.fq.gz ${ meta.trimmed_read1.name }; \
      cp ${ meta.id }_R2_raw_val_2.fq.gz ${ meta.trimmed_read2.name }" : \
      "cp ${ meta.id }_R1_raw_trimmed.fq.gz ${ meta.trimmed_read1.name }"}

    trim_galore -v > versions.txt
    """
}
