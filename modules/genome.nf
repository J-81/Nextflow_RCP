/*
 * Processes related to genome and annotations
 */

process BUILD_STAR {
  conda "${baseDir}/envs/star.yml"
  storeDir ( params.genomeSubsample ?
              "${params.storeDirPath}/${ params.organism }/subsampled/STAR_${ params.ensembl_version }" :
              "${params.storeDirPath}/${ params.organism }/STAR_${ params.ensembl_version }"
            )
  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genomeFasta), path(genomeGtf)
  output:
    path("STAR_REF")
  script:
    """
STAR --runThreadN ${task.cpus} \
--runMode genomeGenerate \
--limitGenomeGenerateRAM ${ task.memory.toBytes() } \
--genomeSAindexNbases 14 \
--genomeDir STAR_REF \
--genomeFastaFiles ${ genomeFasta } \
--sjdbGTFfile ${ genomeGtf } \
--sjdbOverhang 149
    """

}

/*
 * Aligned reads to both the genome and transciptome generated by STAR
 */


process ALIGN_STAR {
  conda "${baseDir}/envs/star.yml"
  label 'maxCPU'
  label 'big_mem'

  input:
    tuple val(sampleID), path(forward_read), path(reverse_read), path(STAR_INDEX_DIR)
  output:
    tuple val(sampleID), \
          path("${ sampleID }Aligned.sortedByCoord.out.bam"), emit: genomeMapping
    tuple val(sampleID), \
          path("${ sampleID }Aligned.toTranscriptome.out.bam"), emit: transcriptomeMapping

  stub:
    """
    touch "${ sampleID }Aligned.sortedByCoord.out.bam" "${ sampleID }Aligned.toTranscriptome.out.bam"
    """

  script:
    """
    STAR --twopassMode Basic \
    --limitBAMsortRAM ${ task.memory.toBytes() } \
    --outFilterType BySJout \
    --outSAMunmapped Within \
    --genomeDir ${ STAR_INDEX_DIR } \
    --outSAMattributes NH HI AS NM MD MC \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMismatchNmax 999 \
    --outFilterMultimapNmax 20 \
    --alignIntronMin 20 \
    --alignSJoverhangMin 8 \
    --alignMatesGapMax 1000000 \
    --alignIntronMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${ task.cpus } \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --outFileNamePrefix ${ sampleID } \
    --readFilesIn ${ forward_read } ${ reverse_read }
    """

}

process BUILD_RSEM {
  conda "${baseDir}/envs/rsem.yml"
  storeDir ( params.genomeSubsample ?
              "${params.storeDirPath}/${ params.organism }/subsampled/RSEM_${ params.ensembl_version }" :
              "${params.storeDirPath}/${ params.organism }/RSEM_${ params.ensembl_version }"
            )

  input:
    tuple path(genomeFasta), path(genomeGtf)
  output:
    path("RSEM_REF")
  script:
    """
    mkdir RSEM_REF
    rsem-prepare-reference --gtf $genomeGtf $genomeFasta RSEM_REF/
    """

}

process COUNT_ALIGNED {
  conda "${baseDir}/envs/rsem.yml"

  input:
    tuple val(sampleID), path(transcriptomeMapping), path(RSEM_REF)
  output:
    tuple val(sampleID), path("${ sampleID }.genes.results"), emit: countsPerGene
    tuple val(sampleID), path("${ sampleID }.isoforms.results"), emit: countsPerIsoform
    tuple val(sampleID), path("${ sampleID }.stat"), emit: stats

  stub:
    """
    touch "${ sampleID }.genes.results" \
          "${ sampleID }.isoforms.results" \
          "${ sampleID }.stat"
    """

  script:
    """
    rsem-calculate-expression --num-threads $task.cpus \
      --paired-end \
      --bam \
      --alignments \
      --no-bam-output \
      --estimate-rspd \
      --seed 12345 \
      --strandedness reverse \
      $transcriptomeMapping \
      $RSEM_REF/ \
      $sampleID
    """

}


/*
 * Download and decompress genome and annotation files
 */

process SUBSAMPLE_GENOME {
  conda "${baseDir}/envs/samtools.yml"
  storeDir "${params.storeDirPath}/ensembl/${params.ensembl_version}/${params.organism}"

  input:
    tuple path(genome_fasta), path(genome_gtf)
  output:
    tuple path("subsampled/${params.genomeSubsample}/${genome_fasta}"), \
          path("subsampled/${params.genomeSubsample}/${genome_gtf}")
  script:
    """
    mkdir -p subsampled/${params.genomeSubsample}
    grep -P "^#|^${params.genomeSubsample}\t" ${genome_gtf} > subsampled/${params.genomeSubsample}/${genome_gtf}

    samtools faidx ${genome_fasta} ${params.genomeSubsample} > subsampled/${params.genomeSubsample}/${genome_fasta}

    """
}

process CONCAT_ERCC {
  storeDir ( params.genomeSubsample ?
              "${params.storeDirPath}/ensembl/${params.ensembl_version}/${params.organism}/ERCC/subsampled" :
              "${params.storeDirPath}/ensembl/${params.ensembl_version}/${params.organism}/ERCC"
              )

  input:
    tuple path(genome_fasta), path(genome_gtf)
    tuple path(ercc_fasta), path(ercc_gtf)
  output:
    tuple path("${params.organism}_and_ERCC.fa"), \
          path("${params.organism}_and_ERCC.gtf")

  script:
  """
  cat ${genome_fasta} ${ercc_fasta} > ${params.organism}_and_ERCC.fa
  cat ${genome_gtf} ${ercc_gtf} > ${params.organism}_and_ERCC.gtf
  """
}
