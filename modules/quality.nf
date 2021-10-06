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
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.PublishTo }",
    mode: params.publish_dir_mode,
    pattern: "*_multiqc_report**"

  label "fastLocal"

  input:
    path("fastqc/*") // any number of fastqc files
  output:
    path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc.html"), emit: html
    path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc_data"), emit: data
    path("versions.txt"), emit: version

  script:
    """
    multiqc -o ${ params.MQCLabel }_multiqc_report -n ${ params.MQCLabel }_multiqc fastqc

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
    echo cutadapt version:\$(cutadapt --version) >> versions.txt
    """
}


process INFER_EXPERIMENT {
  tag "Sample:${ meta.id }"
  // publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.PublishTo }",
  //   mode: params.publish_dir_mode,
  //   pattern: "*_multiqc_report**"

  input:
    tuple val(meta), path(bam_file), path(bed_file) // bam file sorted by coordinate

  output:
    tuple val(meta), path("${ meta.id }_infer_experiment.out"), emit: infer_expt
    path("versions.txt"), emit: version

  script:
    """
    samtools index -@ ${ task.cpus  } ${ bam_file }

    infer_experiment.py -r ${ bed_file } -i ${ bam_file } -s ${ params.quality.rseqc_sample_count } > ${ meta.id }_infer_experiment.out


    infer_experiment.py --version > versions.txt
    """
}


process ASSESS_STRANDEDNESS {
  //tag "Dataset: ${ params.gldsAccession }"
  // publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.PublishTo }",
  //   mode: params.publish_dir_mode,
  //   pattern: "*_multiqc_report**"

  input:
    path("infer_out/*") // a collection of infer_experiment stdout files

  output:
    path("result.txt")

  script:
    """
    assess_strandedness.py infer_out
    """
}
