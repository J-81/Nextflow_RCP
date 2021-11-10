/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  // FastQC performed on reads
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.PublishTo }",
    mode: params.publish_dir_mode,
    pattern: "*_fastqc*"

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

process MULTIQC {
  //tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.PublishTo }",
    mode: params.publish_dir_mode,
    pattern: "*_multiqc**"

  label "fastLocal"

  input:
    path("samples.txt")
    path("mqc_in/*") // any number of multiqc compatible files
    path(multiqc_config)

  output:
    path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc.html"), emit: html
    path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc_data"), emit: data
    path("${ params.MQCLabel }_multiqc_report.zip"), emit: zipped_report
    path("versions.txt"), emit: version

  script:
    config_arg =  multiqc_config.name != "NO_FILE" ? "--config ${ multiqc_config }" : ""
    """
    multiqc --sample-names samples.txt  \
            --interactive -o ${ params.MQCLabel }_multiqc_report \
            -n ${ params.MQCLabel }_multiqc mqc_in \
            ${ config_arg }

    zip -r '${ params.MQCLabel }_multiqc_report.zip' '${ params.MQCLabel }_multiqc_report'

    multiqc --version > versions.txt
    """
}

process TRIMGALORE {
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-TG_Preproc/Fastq",
    mode: params.publish_dir_mode,
    pattern: "*trimmed.fastq.gz"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-TG_Preproc/Trimming_Reports",
    mode: params.publish_dir_mode,
    pattern: "*_trimming_report.txt"

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    path("${ meta.id }*.txt"), emit: reports
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
    echo cutadapt version:\$(cutadapt --version) >> versions.txt
    """
}

