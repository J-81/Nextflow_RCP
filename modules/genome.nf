/*
 * Processes related to genome and annotations
 */

process BUILD_STAR {
  conda "${baseDir}/envs/star.yml"
  tag "Org.:${ meta.organism_sci }  Ensembl.V:${params.ensemblVersion} MaxReadLength:${ max_read_length } GenomeSubsample: ${ params.genomeSubsample }"
  storeDir ( params.genomeSubsample ?
              "${ params.storeDirPath }/${ meta.organism_sci }/readlength_${ max_read_length }/subsampled/${ params.genomeSubsample }/STAR_ensembl_${ params.ensemblVersion }" :
              "${ params.storeDirPath }/${ meta.organism_sci }/readlength_${ max_read_length }/STAR_ensembl_${ params.ensemblVersion }"
            )
  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)
  output:
    path("STAR_REF"), emit: build
    path("versions.txt"), emit: version

  script:
    max_read_length = meta.paired_end ? [meta.read_length_R1, meta.read_length_R2].max() : meta.read_length_R1
    if (!max_read_length) { throw new Exception("NullOrFalse Max Read Length: ${max_read_length}") }
    """
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM ${ task.memory.toBytes() } \
    --genomeSAindexNbases 14 \
    --genomeDir STAR_REF \
    --genomeFastaFiles ${ genomeFasta } \
    --sjdbGTFfile ${ genomeGtf } \
    --sjdbOverhang ${ max_read_length - 1 }

    echo STAR_version: `STAR --version` > versions.txt
    """

}

/*
 * Aligned reads to both the genome and transciptome generated by STAR
 */


process ALIGN_STAR {
  conda "${baseDir}/envs/star.yml"
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }"
  label 'maxCPU'
  label 'big_mem'

  input:
    tuple val( meta ), path( reads ), path(STAR_INDEX_DIR)
  output:
    tuple val(meta), path("${ meta.STAR_Alignment_dir }")

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
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --runThreadN ${ task.cpus } \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --outFileNamePrefix '${ meta.STAR_Alignment_dir }/${ meta.id }_' \
    --readFilesIn ${ reads }
    """

}

process BUILD_RSEM {
  conda "${baseDir}/envs/rsem.yml"
  tag "Organism: ${ meta.organism_sci }  Ensembl Version: ${params.ensemblVersion}"
  storeDir ( params.genomeSubsample ?
              "${params.storeDirPath}/${ meta.organism_sci }/subsampled/${ params.genomeSubsample }/RSEM_ensembl_${ params.ensemblVersion }" :
              "${params.storeDirPath}/${ meta.organism_sci }/RSEM_ensembl_${ params.ensemblVersion }"
            )

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)
  output:
    path("RSEM_REF"), emit: build
    path("versions.txt"), emit: version


  script:
    """
    mkdir RSEM_REF
    rsem-prepare-reference --gtf $genomeGtf $genomeFasta RSEM_REF/

    rsem-calculate-expression --version > versions.txt
    """

}

process COUNT_ALIGNED {
  conda "${baseDir}/envs/rsem.yml"
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }"

  input:
    tuple val(meta), path("starOutput/*"), path(RSEM_REF)
  output:

    tuple val(meta), path("${ meta.RSEM_Counts_dir }/*")

  script:
    """
    rsem-calculate-expression --num-threads $task.cpus \
      ${ meta.paired_end ? '--paired-end' : '' } \
      --bam \
      --alignments \
      --no-bam-output \
      --estimate-rspd \
      --seed 12345 \
      --strandedness reverse \
      starOutput/${meta.id}/${meta.id}_Aligned.toTranscriptome.out.bam \
      $RSEM_REF/ \
      $meta.id

    # move into output directory
    mkdir tmp
    mv ${meta.id}* tmp
    mkdir -p ${ meta.RSEM_Counts_dir }
    mv tmp/* ${ meta.RSEM_Counts_dir }
    """
}

process QUANTIFY_GENES {
  conda "${baseDir}/envs/RNAseq_Rtools.yml"
  tag "Dataset: ${ params.gldsAccession }"

  input:
    path("samples.txt")
    path("03-RSEM_COUNTS/*.genes.results")

  output:
    tuple path("RSEM_Unnormalized_Counts.csv"), path("NumNonZeroGenes.csv")

  script:
    """
    Quantitate_non-zero_genes_per_sample.R
    """

}


/*
 * Download and decompress genome and annotation files
 */

process SUBSAMPLE_GENOME {
  conda "${baseDir}/envs/samtools.yml"
  storeDir "${params.storeDirPath}/ensembl/${params.ensemblVersion}/${ meta.organism_sci }"

  input:
    tuple path(genome_fasta), path(genome_gtf)
    val(meta)
  output:
    tuple path("subsampled/${params.genomeSubsample}/${genome_fasta}"), \
          path("subsampled/${params.genomeSubsample}/${genome_gtf}"), emit: build
    path("versions.txt"), emit: version

  script:
    """
    mkdir -p subsampled/${params.genomeSubsample}
    grep -P "^#|^${params.genomeSubsample}\t" ${genome_gtf} > subsampled/${params.genomeSubsample}/${genome_gtf}

    samtools faidx ${genome_fasta} ${params.genomeSubsample} > subsampled/${params.genomeSubsample}/${genome_fasta}

    samtools --version > versions.txt
    """
}

process CONCAT_ERCC {
  errorStrategy 'retry'
  maxRetries 3
  storeDir ( params.genomeSubsample ?
              "${params.storeDirPath}/ensembl/${params.ensemblVersion}/${ meta.organism_sci }/ERCC/subsampled" :
              "${params.storeDirPath}/ensembl/${params.ensemblVersion}/${ meta.organism_sci }/ERCC"
              )

  input:
    tuple path(genome_fasta), path(genome_gtf)
    tuple path(ercc_fasta), path(ercc_gtf)
    val(meta)
  output:
    tuple path("${ meta.organism_sci }_and_ERCC.fa"), \
          path("${ meta.organism_sci }_and_ERCC.gtf")

  when:
    meta.has_ercc

  script:
  """
  cat ${genome_fasta} ${ercc_fasta} > ${ meta.organism_sci }_and_ERCC.fa
  cat ${genome_gtf} ${ercc_gtf} > ${ meta.organism_sci }_and_ERCC.gtf
  """
}
