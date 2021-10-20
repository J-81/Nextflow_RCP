// Workflow that determines the strandedness of reads compared to a reference genome bed file
include { RSEQC_ALL } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/logs")
include { ASSESS_STRANDEDNESS } from './modules/rseqc.nf'


workflow strandedness{
  take:
    bam_array // array: sample-wise tuples (meta, bam_file)
    genome_bed 
 
  main:
     
     bam_array | combine( genome_bed ) |  set { ch_bam_bed }
     RSEQC_ALL( ch_bam_bed ) 

     RSEQC_ALL.out.infer_expt | map { it[1] }
                                     | collect
                                     | set { ch_infer_expt }
    
     ch_infer_expt | ASSESS_STRANDEDNESS
     RSEQC_ALL.out.all | map{ it[1..it.size()-1] } 
                       | collect 
                       | set { ch_rseqc_logs }
     

  emit:
     strandedness = ASSESS_STRANDEDNESS.out 
     infer_expt = ch_infer_expt
     versions = RSEQC_ALL.out.version
     rseqc_logs = ch_rseqc_logs
}
