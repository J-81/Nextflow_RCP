/*
 * Processes related to genome and annotations
 */

process BUILD_STAR {
  // Builds STAR index, this is ercc-spike-in, organism, read length and ensembl version specific
  tag "Refs:${ genomeFasta },${ genomeGtf }, Ensembl.V:${params.ensemblVersion} MaxReadLength:${ max_read_length } GenomeSubsample: ${ params.genomeSubsample }"
  storeDir "${ params.derivedStorePath }/STAR_Indices/${ params.ref_source }_release${params.ensemblVersion}/${ meta.organism_sci.capitalize() }"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)
    val(max_read_length) // Based on fastQC report for all samples

  output:
    path("${ genomeFasta.baseName }_RL-${ max_read_length.toInteger() }"), emit: build
    path("${ genomeFasta.baseName }_RL-${ max_read_length.toInteger() }/genomeParameters.txt") // Check for completion, only successful builds should generate this file, this is required as the process error is NOT currently used to raised an exception in the python wrapper.

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
--genomeDir ${ genomeFasta.baseName }_RL-${ max_read_length.toInteger() } \
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
    --genomeDir ${ genomeFasta.baseName }_RL-${ max_read_length.toInteger() } \
    --genomeFastaFiles ${ genomeFasta } \
    --sjdbGTFfile ${ genomeGtf } \
    --sjdbOverhang ${ max_read_length.toInteger() - 1 }
    '''
    print(command)
    command = command.format(suggested=suggested)
    
    output_rerun = subprocess.check_output(shlex.split(command), shell=False)
    print(output_rerun)

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
  // TODO: make '--alignMatesGapMax 1000000' conditional on PE
  tag "Sample: ${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/02-STAR_Alignment",
    mode: params.publish_dir_mode,
    pattern: "${ meta.id }/${ meta.id }*"
  label 'maxCPU'
  label 'big_mem'

  input:
    tuple val( meta ), path( reads ), path(STAR_INDEX_DIR)

  output:
    path("${ meta.id }/${ meta.id }*"), emit: publishables // used to ensure direct files are available for publishing directive
    path("${ meta.id }/${ meta.id}_Log.final.out"), emit: alignment_logs
    tuple val(meta), path("${ meta.id }/${ meta.id }_Aligned.sortedByCoord.out.bam"), emit: bam_by_coord
    tuple val(meta), path("${ meta.id }/${ meta.id }_Aligned.toTranscriptome.out.bam"), emit: bam_to_transcriptome
    path("${ meta.id }/${ meta.id }_ReadsPerGene.out.tab"), emit: read_per_gene
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
    --quantMode TranscriptomeSAM GeneCounts\
    --outFileNamePrefix '${ meta.id }/${ meta.id }_' \
    --readFilesIn ${ reads }

    echo ALIGN_STAR_version: `STAR --version` > versions.txt
    """

}

process BUILD_RSEM {
  // Builds RSEM index, this is ercc-spike-in, organism, and ensembl version specific
  tag "Refs:${ genomeFasta },${ genomeGtf }, Ensembl Version: ${params.ensemblVersion}, GenomeSubsample: ${ params.genomeSubsample }"
  storeDir "${ params.derivedStorePath }/RSEM_Indices/${ params.ref_source }_release${params.ensemblVersion}/${ meta.organism_sci.capitalize() }"

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)

  output:
    path("${ genomeFasta.baseName }"), emit: build
    path("${ genomeFasta.baseName }/${ organism_str }.grp") // to ensure check expected file contents exist


  script:
    // e.g. 'mus_musculus' should become 'Mmus'
    organism_str = "${ meta.organism_sci.substring(0,1).toUpperCase() }${ meta.organism_sci.split('_')[1].substring(0,3).toLowerCase() }"
    """
    mkdir  ${ genomeFasta.baseName }
    rsem-prepare-reference --gtf $genomeGtf $genomeFasta ${ genomeFasta.baseName }/${ organism_str }

    # echo Build_RSEM_version: `rsem-calculate-expression --version` > versions.txt
    """

}

process COUNT_ALIGNED {
  // Generates gene and isoform counts from alignments
  tag "Sample: ${ meta.id }, strandedness: ${ strandedness } "
  publishDir "${ params.outputDir }/${ params.gldsAccession }/03-RSEM_Counts",
    mode: params.publish_dir_mode,
    pattern: "${ meta.id }*"

  input:
    tuple val(meta), path("${meta.id}_Aligned.toTranscriptome.out.bam"), path(RSEM_REF)
    val(strandedness)

  output:
    tuple val(meta), path("${ meta.id }*"), emit: counts
    path("versions.txt"), emit: version

  script:
    strandedness_opt_map = ["sense":"forward","antisense":"reverse","unstranded":"none"]
    // e.g. 'mus_musculus' should become 'Mmus'
    organism_str = "${ meta.organism_sci.substring(0,1).toUpperCase() }${ meta.organism_sci.split('_')[1].substring(0,3).toLowerCase() }"
    """
    rsem-calculate-expression --num-threads $task.cpus \
      ${ meta.paired_end ? '--paired-end' : '' } \
      --bam \
      --alignments \
      --no-bam-output \
      --estimate-rspd \
	    --seed-length 20 \
      --seed 12345 \
      --strandedness ${ strandedness_opt_map.get(strandedness) } \
      ${meta.id}_Aligned.toTranscriptome.out.bam \
      ${ RSEM_REF }/${ organism_str } \
      ${ meta.id }

    echo COUNT_RSEM_version: `rsem-calculate-expression --version` > versions.txt
    """
}

process QUANTIFY_RSEM_GENES {
  // An R script that extracts gene counts by sample to a table
  publishDir "${ params.outputDir }/${ params.gldsAccession }/03-RSEM_Counts",
    mode: params.publish_dir_mode

  input:
    path("samples.txt")
    path("03-RSEM_Counts/*")

  output:
    tuple path("RSEM_Unnormalized_Counts.csv"), path("RSEM_NumNonZeroGenes.csv")

  script:
    """
    Quantitate_non-zero_genes_per_sample.R
    """

}

process QUANTIFY_STAR_GENES {
  // An R script that extracts gene counts by sample to a table
  publishDir "${ params.outputDir }/${ params.gldsAccession }/02-STAR_Alignment",
    mode: params.publish_dir_mode

  input:
    path("samples.txt")
    path("02-STAR_Alignment/*")
    val(strandedness)

  output:
    tuple path("STAR_Unnormalized_Counts.csv"), path("STAR_NumNonZeroGenes.csv")

  script:
    """
    Quantitate_non-zero_genes_per_sample_STAR.R ${strandedness}
    """

}

process SUBSAMPLE_GENOME {
  // Extracts a user-specified sequence from the larger reference fasta and gtf file
  tag "Sequence:'${ params.genomeSubsample }'"
  storeDir "${params.derivedStorePath}/subsampled_files/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"

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
  storeDir "${params.referenceStorePath}/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
          

  input:
    tuple path(genome_fasta), path(genome_gtf)
    tuple path(ercc_fasta), path(ercc_gtf)
    val(organism_sci)
    val(has_ercc)

  output:
    tuple path("${ genome_fasta.baseName }_and_ERCC92.fa"), \
          path("${ genome_gtf.baseName }_and_ERCC92.gtf")

  when:
    has_ercc

  script:
  """
  cat ${genome_fasta} ${ercc_fasta} > ${ genome_fasta.baseName }_and_ERCC92.fa
  cat ${genome_gtf} ${ercc_gtf} > ${ genome_gtf.baseName }_and_ERCC92.gtf
  """
}

process TO_PRED {
  // Converts reference gtf into pred 
  storeDir "${ params.derivedStorePath }/Genome_GTF_BED_Files/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
          

  input:
    path(genome_gtf)
    val(organism_sci)

  output:
    path("${ genome_gtf }.genePred")

  script:
  """
  gtfToGenePred -geneNameAsName2 ${ genome_gtf } ${ genome_gtf }.genePred
  """
}


process TO_BED {
  // Converts reference genePred into Bed format
  storeDir "${ params.derivedStorePath }/Genome_GTF_BED_Files/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
          

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
