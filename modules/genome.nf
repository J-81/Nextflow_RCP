/*
 * Processes related to genome and annotations
 */

process BUILD_STAR {
  // Builds STAR index, this is ercc-spike-in, organism, read length and ensembl version specific
  tag "Refs:${ genomeFasta },${ genomeGtf }, Ensembl.V:${params.ensemblVersion} MaxReadLength:${ max_read_length } GenomeSubsample: ${ params.genomeSubsample }"
  storeDir "${ params.storeDirPath }/STAR_Indices"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)
    val(max_read_length) // Based on fastQC report for all samples

  output:
    path("STAR_REF_${ genomeFasta.baseName }"), emit: build
    path("STAR_REF_${ genomeFasta.baseName }/Log.out") // Check for log completed

  script:
    """
#! /usr/bin/env python
    
import subprocess
import shlex

command = '''
STAR --runThreadN ${task.cpus} \
--runMode genomeGenerate \
--limitGenomeGenerateRAM ${ task.memory.toBytes() } \
--genomeSAindexNbases 14 \
--genomeDir STAR_REF_${ genomeFasta.baseName } \
--genomeFastaFiles ${ genomeFasta } \
--sjdbGTFfile ${ genomeGtf } \
--sjdbOverhang ${ max_read_length.toInteger() - 1 }
'''

def rerun(suggested):
    command = '''
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM ${ task.memory.toBytes() } \
    --genomeSAindexNbases {suggested} \
    --genomeDir STAR_REF_${ genomeFasta.baseName } \
    --genomeFastaFiles ${ genomeFasta } \
    --sjdbGTFfile ${ genomeGtf } \
    --sjdbOverhang ${ max_read_length.toInteger() - 1 }
    '''
    print(command)
    command = command.format(suggested=suggested)
    
    subprocess.Popen(shlex.split(command), shell=False)

# invoke initial build process
process = subprocess.Popen(shlex.split(command), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# Poll process.stdout
while True:
    output = process.stderr.readline()
    if process.poll() is not None:
        break
    if output:
        line = output.strip().decode()
        print(line)
        if "!!!!! WARNING: --genomeSAindexNbases" in line: # indicating better parameters
            suggested = int(line.split()[-1])
            print(f"Restarting build with suggested '--genomeSAindexNbases' parameter ({suggested})")
            process.kill()
            rerun(suggested)
            break
        
    rc = process.poll()
    """

}


process ALIGN_STAR {
  // Aligns reads against STAR index
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    mode: params.publish_dir_mode,
    pattern: "${ meta.STAR_Alignment_dir }"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple val( meta ), path( reads ), path(STAR_INDEX_DIR)

  output:
    tuple val(meta), path("${ meta.STAR_Alignment_dir }"), emit: alignments
    path("${ meta.STAR_Alignment_dir }/${ meta.id}_Log.final.out"), emit: alignment_logs
    tuple val(meta), path("${ meta.STAR_Alignment_dir }/${ meta.id }_Aligned.sortedByCoord.out.bam"), emit: bam_by_coord
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
  tag "Refs:${ genomeFasta },${ genomeGtf }, Ensembl Version: ${params.ensemblVersion}, GenomeSubsample: ${ params.genomeSubsample }"
  storeDir "${ params.storeDirPath }/RSEM_Indices"

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)

  output:
    path("RSEM_REF_${ genomeFasta.baseName }"), emit: build
    path("RSEM_REF_${ genomeFasta.baseName }/.grp") // to ensure check expected file contents exist


  script:
    """
    mkdir  RSEM_REF_${ genomeFasta.baseName }
    rsem-prepare-reference --gtf $genomeGtf $genomeFasta RSEM_REF_${ genomeFasta.baseName }/

    # echo Build_RSEM_version: `rsem-calculate-expression --version` > versions.txt
    """

}

process COUNT_ALIGNED {
  // Generates gene and isoform counts from alignments
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    mode: params.publish_dir_mode,
    pattern: "${ meta.RSEM_Counts_dir }/*"

  input:
    tuple val(meta), path("starOutput/*"), path(RSEM_REF)
    val(strandedness)

  output:
    tuple val(meta), path("${ meta.RSEM_Counts_dir }/*"), emit: counts
    path("versions.txt"), emit: version

  script:
    strandedness_opt_map = ["sense":"forward","antisense":"reverse","unstranded":"none"]
    """
    rsem-calculate-expression --num-threads $task.cpus \
      ${ meta.paired_end ? '--paired-end' : '' } \
      --bam \
      --alignments \
      --no-bam-output \
      --estimate-rspd \
      --seed 12345 \
      --strandedness ${ strandedness_opt_map.get(strandedness) } \
      starOutput/${meta.id}/${meta.id}_Aligned.toTranscriptome.out.bam \
      ${ RSEM_REF }/ \
      ${ meta.id }

    # move results into output directory
    mkdir tmp
    mv ${ meta.id }* tmp
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
  tag "Sequence:'${ params.genomeSubsample }'"
  storeDir "${params.storeDirPath}/ensembl/version-${params.ensemblVersion}/${ organism_sci }"

  input:
    tuple path(genome_fasta), path(genome_gtf)
    val(organism_sci)

  output:
    tuple path("${ genome_fasta.baseName }_sub_${ params.genomeSubsample  }.fa"), \
          path("${ genome_gtf.baseName }_sub_${ params.genomeSubsample }.gtf"), emit: build

  script:
    """
    samtools faidx ${genome_fasta} ${params.genomeSubsample} > ${ genome_fasta.baseName }_sub_${ params.genomeSubsample }.fa

    grep -P "^#|^${params.genomeSubsample}\t" ${genome_gtf} > ${ genome_gtf.baseName }_sub_${ params.genomeSubsample  }.gtf
    """
}

process CONCAT_ERCC {
  // Concanates ERCC fasta and gtf to reference fasta and gtf
  errorStrategy 'retry'
  maxRetries 3 // This addresses a very rare unexpected error where the command finishes but output is not produced.
  storeDir "${params.storeDirPath}/ensembl/version-${params.ensemblVersion}/${ organism_sci }"
          

  input:
    tuple path(genome_fasta), path(genome_gtf)
    tuple path(ercc_fasta), path(ercc_gtf)
    val(organism_sci)
    val(has_ercc)

  output:
    tuple path("${ genome_fasta.baseName }_and_ERCC.fa"), \
          path("${ genome_gtf.baseName }_and_ERCC.gtf")

  when:
    has_ercc

  script:
  """
  cat ${genome_fasta} ${ercc_fasta} > ${ genome_fasta.baseName }_and_ERCC.fa
  cat ${genome_gtf} ${ercc_gtf} > ${ genome_gtf.baseName }_and_ERCC.gtf
  """
}

process TO_PRED {
  // Converts reference gtf into pred 
          

  input:
    path(genome_gtf)

  output:
    path("${ genome_gtf }.genePred")

  script:
  """
  gtfToGenePred -geneNameAsName2 ${ genome_gtf } ${ genome_gtf }.genePred
  """
}


process TO_BED {
  // Converts reference genePred into Bed format
  storeDir "${params.storeDirPath}/ensembl/version-${params.ensemblVersion}/${ organism_sci }"
          

  input:
    path(genome_pred)
    val(organism_sci)

  output:
    path("${ genome_pred.baseName }.bed")

  script:
  """
  genePredToBed ${ genome_pred } ${ genome_pred.baseName }.bed
  """
}
