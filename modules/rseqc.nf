process ASSESS_STRANDEDNESS {

  input:
    path("infer_out/*") // a collection of infer_experiment stdout files

  output:
    path("result.txt")

  script:
    """
    assess_strandedness.py infer_out
    """
}

process RSEQC_ALL {
  tag "Sample:${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ params.PublishTo }",
     mode: params.publish_dir_mode,
     pattern: "${ meta.id }/*"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bed_file) // bam file sorted by coordinate

  output:
    tuple val(meta), path("${ meta.id }/${ meta.id }.infer_experiment_out"), emit: infer_expt
    path("versions.txt"), emit: version
    tuple val(meta), path("${ meta.id }/${ meta.id }.infer_experiment_out"), \
                     path("${ meta.id }/${ meta.id }.geneBodyCoverage.txt"), \
                     path("${ meta.id }/${ meta.id }.inner_distance_freq.txt"), \
                     path("${ meta.id }/${ meta.id }.read_distribution_out"), \
                     emit: all 

  script:
    sorted_bam_fname = bam_file.name.replaceAll('.out.bam','_sorted.out.bam')
    mem_MB_per_thread = task.memory.toMega().intdiv(task.cpus)
    """    
    mkdir ${ meta.id }
    samtools sort -m ${ mem_MB_per_thread  }M \
                  --threads ${ task.cpus } \
                  -o ${ sorted_bam_fname } \
                  ${ bam_file }

    samtools index -@ ${ task.cpus  } ${ sorted_bam_fname }

    infer_experiment.py -r ${ bed_file } -i ${ sorted_bam_fname } -s ${ params.quality.rseqc_sample_count } > ${ meta.id }/${ meta.id }.infer_experiment_out

    geneBody_coverage.py -r ${ bed_file} -i ${ sorted_bam_fname } -o ${ meta.id }/${ meta.id }

    inner_distance.py -r ${ bed_file } -i ${ sorted_bam_fname } -k 15000000 -l -150 -u 350 -o ${ meta.id }/${ meta.id } // run for single end is nonsense and generates blank file, okay for now, make sure to disable output to versions
    read_distribution.py -r ${ bed_file } -i ${ sorted_bam_fname } > ${ meta.id }/${ meta.id }.read_distribution_out 

    
    # VERSIONS
    echo "RSeQC tools versions below:\n" > versions.txt 
    infer_experiment.py --version >> versions.txt
    geneBody_coverage.py --version >> versions.txt
    ${ meta.paired_end ? 'inner_distance.py --version >> versions.txt' : '' }
    read_distribution.py  --version >> versions.txt
    """
}

