/*
 * Processes related to genome and annotations
 */

process BUILD_STAR {
  // Builds STAR index, this is ercc-spike-in, organism, read length and ensembl version specific
  storeDir "${ params.derivedStorePath }/STAR_Indices/${ params.ref_source }_release${params.ensemblVersion}/${ meta.organism_sci.capitalize() }"
  tag "storeDir: ...${ task.storeDir.toString()[-params.shorten_storeDir_tag..-1] } Target(s): ${ build_dir }, ${ build_dir }/genomeParameters.txt"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)
    val(max_read_length) // Based on fastQC report for all samples

  output:
    path("${ build_dir }"), emit: build
    path("${ build_dir }/genomeParameters.txt") // Check for completion, only successful builds should generate this file, this is required as the process error is NOT currently used to raised an exception in the python wrapper.

  script:
    build_dir = "${ genomeFasta.baseName.split('\\.')[0] }_RL-${ max_read_length.toInteger() }"
    """
#! /usr/bin/env python
    
import subprocess
import shlex

command = '''
STAR --runThreadN ${task.cpus} \
--runMode genomeGenerate \
--limitGenomeGenerateRAM ${ task.memory.toBytes() } \
--genomeSAindexNbases 14 \
--genomeDir ${ build_dir } \
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
    --genomeDir ${ build_dir } \
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
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    mode: params.publish_dir_mode,
    pattern: "${ meta.STAR_Alignment_dir }/${ meta.id }**"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple val( meta ), path( reads ), path(STAR_INDEX_DIR)

  output:
    tuple val(meta), path("${ meta.STAR_Alignment_dir }"), emit: alignments
    path("${ meta.STAR_Alignment_dir }/${ meta.id}_Log.final.out"), emit: alignment_logs
    tuple val(meta), path("${ meta.STAR_Alignment_dir }/${ meta.id }_Aligned.sortedByCoord.out.bam"), emit: bam_by_coord
    path("${ meta.STAR_Alignment_dir }/${ meta.id }**")
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
  storeDir "${ params.derivedStorePath }/RSEM_Indices/${ params.ref_source }_release${params.ensemblVersion}"
  tag "storeDir: ...${ task.storeDir.toString()[-params.shorten_storeDir_tag..-1] } Target(s): ${ build_prefix }*"

  input:
    tuple path(genomeFasta), path(genomeGtf)
    val(meta)

  output:
    path("${ build_dir }"), emit: build
    path("${ build_prefix }.grp") // to ensure check expected file contents exist


  script:
    ercc_substring = meta.has_ercc ? '_w_ERCC' : ''
    build_dir = "${ meta.organism_sci.capitalize() }${ ercc_substring }"
    organism_substring = "${ meta.organism_sci.capitalize()[0] }${ meta.organism_sci.split('_')[1][0..2] }${ ercc_substring }"
    build_prefix = "${ build_dir }/${ organism_substring }"
    """
    mkdir  ${ build_dir }
    rsem-prepare-reference --gtf $genomeGtf $genomeFasta $build_prefix 

    chmod g-w ${ build_prefix }*

    # echo Build_RSEM_version: `rsem-calculate-expression --version` > versions.txt
    """

}

process COUNT_ALIGNED {
  // Generates gene and isoform counts from alignments
  tag "Sample: ${ meta.id }, strandedness: ${ strandedness } "
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

    ercc_substring = meta.has_ercc ? '_w_ERCC' : ''
    build_dir = "${ meta.organism_sci.capitalize() }${ ercc_substring }"
    organism_substring = "${ meta.organism_sci.capitalize()[0] }${ meta.organism_sci.split('_')[1][0..2] }${ ercc_substring }"
    build_prefix = "${ build_dir }/${ organism_substring }"
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
      ${ build_prefix } \
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
  storeDir "${params.derivedStorePath}/subsampled_files/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
  tag "storeDir: ...${ task.storeDir.toString()[-params.shorten_storeDir_tag..-1] } Target(s): ${ genome_fasta.baseName }_sub_${ params.genomeSubsample  }.fa, ${ genome_gtf.baseName }_sub_${ params.genomeSubsample }.gtf"

  input:
    tuple path(genome_fasta), path(genome_gtf)
    val(organism_sci)

  output:
    tuple path("${ genome_fasta.baseName }_sub_${ params.genomeSubsample  }.fa"), \
          path("${ genome_gtf.baseName }_sub_${ params.genomeSubsample }.gtf"), emit: build

  script:
    """
    # ensure proper permissions for generated file
    umask 0022
  
    samtools faidx ${genome_fasta} ${params.genomeSubsample} > ${ genome_fasta.baseName }_sub_${ params.genomeSubsample }.fa

    grep -P "^#|^${params.genomeSubsample}\t" ${genome_gtf} > ${ genome_gtf.baseName }_sub_${ params.genomeSubsample  }.gtf
    """
}

process CONCAT_ERCC {
  // Concanates ERCC fasta and gtf to reference fasta and gtf
  errorStrategy 'retry'
  maxRetries 3 // This addresses a very rare unexpected error where the command finishes but output is not produced.
  storeDir "${params.referenceStorePath}/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
  tag "storeDir: ...${ task.storeDir.toString()[-params.shorten_storeDir_tag..-1] } Target(s): ${ genome_fasta.baseName }_and_${ ercc_fasta.name }, ${ genome_gtf.baseName }_and_${ ercc_gtf.name }"
          

  input:
    tuple path(genome_fasta), path(genome_gtf)
    tuple path(ercc_fasta), path(ercc_gtf)
    val(organism_sci)
    val(has_ercc)

  output:
    tuple path("${ genome_fasta.baseName }_and_${ ercc_fasta.name }"), \
          path("${ genome_gtf.baseName }_and_${ ercc_gtf.name }")

  when:
    has_ercc

  script:
  """
  # ensure proper permissions for generated file
  umask 0022
  
  cat ${genome_fasta} ${ercc_fasta} > ${ genome_fasta.baseName }_and_${ ercc_fasta.name }
  cat ${genome_gtf} ${ercc_gtf} > ${ genome_gtf.baseName }_and_${ ercc_gtf.name }
  """
}

process TO_PRED {
  // Converts reference gtf into pred 
  tag "storeDir Target: ...${ task.storeDir.toString()[-params.shorten_storeDir_tag..-1] }/${ genome_gtf }.genePred"
  storeDir "${ params.derivedStorePath }/Genome_GTF_BED_Files/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
          

  input:
    path(genome_gtf)
    val(organism_sci)

  output:
    path("${ genome_gtf }.genePred")

  script:
  """
  # ensure proper permissions for generated file
  umask 0022
  
  gtfToGenePred -geneNameAsName2 ${ genome_gtf } ${ genome_gtf }.genePred
  """
}


process TO_BED {
  // Converts reference genePred into Bed format
  tag "storeDir Target: ...${ task.storeDir.toString()[-params.shorten_storeDir_tag..-1] }/${ genome_pred.baseName }.bed"
  storeDir "${ params.derivedStorePath }/Genome_GTF_BED_Files/${ params.ref_source }_release${params.ensemblVersion}/${ organism_sci.capitalize() }"
          

  input:
    path(genome_pred)
    val(organism_sci)

  output:
    path("${ genome_pred.baseName }.bed")

  script:
  """
  # ensure proper permissions for generated file
  umask 0022
  
  genePredToBed ${ genome_pred } ${ genome_pred.baseName }.bed
  """
}
