/*
 * Processes related to genome and annotations
 */

// A string used to modify the filename for subsampled indices
subsample_mod = params.genomeSubsample ? "_subsampled-${ params.genomeSubsample }" : ""

process BUILD_STAR {
  // Builds STAR index, this is ercc-spike-in, organism, read length and ensembl version specific
  tag "Org.:${ meta.organism_sci }  Ensembl.V:${params.ensemblVersion} MaxReadLength:${ max_read_length } GenomeSubsample: ${ params.genomeSubsample }"
  storeDir "${ params.storeDirPath }/STAR_Indices/ensembl_release${params.ensemblVersion}/${ meta.organism_sci }${ ercc_mod }_RL-${ max_read_length }${ subsample_mod }"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)
    val(max_read_length) // Based on fastQC report for all samples

  output:
    path("STAR_REF"), emit: build
    // path("versions.txt"), emit: version

  script:
    // Filename modifier string for indices with ERCC included
    ercc_mod = meta.has_ercc ? "_w_ERCC" : ""
    """
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM ${ task.memory.toBytes() } \
    --genomeSAindexNbases 14 \
    --genomeDir STAR_REF \
    --genomeFastaFiles ${ genomeFasta } \
    --sjdbGTFfile ${ genomeGtf } \
    --sjdbOverhang ${ max_read_length.toInteger() - 1 }

    # echo Build_STAR_version: `STAR --version` > versions.txt
    """

}


process ALIGN_STAR {
  // Aligns reads against STAR index
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    mode: params.publish_dir_mode

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple val( meta ), path( reads ), path(STAR_INDEX_DIR)

  output:
    tuple val(meta), path("${ meta.STAR_Alignment_dir }"), emit: alignments
    path("versions.txt"), emit: version

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

    echo ALIGN_STAR_version: `STAR --version` > versions.txt
    """

}

process BUILD_RSEM {
  // Builds RSEM index, this is ercc-spike-in, organism, and ensembl version specific
  tag "Organism: ${ meta.organism_sci }  Ensembl Version: ${params.ensemblVersion}"
  storeDir "${ params.storeDirPath }/RSEM_Indices/ensembl_release${params.ensemblVersion}/${ meta.organism_sci }${ ercc_mod }${ subsample_mod }"

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)

  output:
    path("RSEM_REF"), emit: build
    // path("versions.txt"), emit: version


  script:
    // Filename modifier string for indices with ERCC included
    ercc_mod = meta.has_ercc ? "_w_ERCC" : ""
    """
    mkdir RSEM_REF
    rsem-prepare-reference --gtf $genomeGtf $genomeFasta RSEM_REF/${meta.organism_sci}${ercc_mod}

    # echo Build_RSEM_version: `rsem-calculate-expression --version` > versions.txt
    """

}

process COUNT_ALIGNED {
  // Generates gene and isoform counts from alignments
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    mode: params.publish_dir_mode

  input:
    tuple val(meta), path("starOutput/*"), path(RSEM_REF)

  output:
    tuple val(meta), path("${ meta.RSEM_Counts_dir }/*"), emit: counts
    path("versions.txt"), emit: version

  script:
    ercc_mod = meta.has_ercc ? "_w_ERCC" : ""
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
      $RSEM_REF/${meta.organism_sci}${ercc_mod} \
      $meta.id

    # move into output directory
    mkdir tmp
    mv ${meta.id}* tmp
    mkdir -p ${ meta.RSEM_Counts_dir }
    mv tmp/* ${ meta.RSEM_Counts_dir }

    echo COUNT_RSEM_version: `rsem-calculate-expression --version` > versions.txt
    """
}

process QUANTIFY_GENES {
  // An R script that extracts gene counts by sample to a table
  tag "Dataset: ${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/03-RSEM_Counts",
    mode: params.publish_dir_mode

  input:
    path("samples.txt")
    path("03-RSEM_Counts/*")

  output:
    tuple path("RSEM_Unnormalized_Counts.csv"), path("NumNonZeroGenes.csv")

  script:
    """
    Quantitate_non-zero_genes_per_sample.R
    """

}

process SUBSAMPLE_GENOME {
  // Extracts a user-specified sequence from the larger reference fasta and gtf file
  storeDir "${params.storeDirPath}/ensembl/${params.ensemblVersion}/${ meta.organism_sci }"

  input:
    tuple path(genome_fasta), path(genome_gtf)
    val(meta)
  output:
    tuple path("subsampled/${params.genomeSubsample}/${genome_fasta}"), \
          path("subsampled/${params.genomeSubsample}/${genome_gtf}"), emit: build

  script:
    """
    mkdir -p subsampled/${params.genomeSubsample}
    grep -P "^#|^${params.genomeSubsample}\t" ${genome_gtf} > subsampled/${params.genomeSubsample}/${genome_gtf}

    samtools faidx ${genome_fasta} ${params.genomeSubsample} > subsampled/${params.genomeSubsample}/${genome_fasta}

    """
}

process CONCAT_ERCC {
  // Concanates ERCC fasta and gtf to reference fasta and gtf
  errorStrategy 'retry'
  maxRetries 3 // This addresses a very rare unexpected error where the command finishes but output is not produced.
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
