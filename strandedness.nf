// Workflow that determines the strandedness of reads compared to a reference genome bed file
include { INFER_EXPERIMENT;
          ASSESS_STRANDEDNESS } from './modules/quality.nf'

include { MULTIQC as INFER_EXPT_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "01-TG_Preproc/RSeQC_Reports", MQCLabel:"infer_expt")

workflow strandedness{
  take:
    bam_array // array: sample-wise tuples (meta, bam_file)
    genome_bed 
 
  main:
     
     bam_array | combine( genome_bed ) |  set { ch_bam_bed }
     INFER_EXPERIMENT( ch_bam_bed ) 

     INFER_EXPERIMENT.out.infer_expt | map { it[1] }
                                     | collect
                                     | set { ch_infer_expt }
     
     ch_infer_expt | (ASSESS_STRANDEDNESS & INFER_EXPT_MULTIQC ) 

  emit:
     strandedness = ASSESS_STRANDEDNESS.out 
     versions = INFER_EXPERIMENT.out.version
}
