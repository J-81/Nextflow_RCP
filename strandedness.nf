// Workflow that determines the strandedness of reads compared to a reference genome bed file
include { SORT_INDEX_BAM } from './modules/rseqc.nf' addParams(PublishTo: "02-STAR_Alignment")
include { INFER_EXPERIMENT } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/infer_experiment_logs")
include { GENEBODY_COVERAGE } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/genebody_coverage_logs")
include { INNER_DISTANCE } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/inner_distance_logs")
include { READ_DISTRIBUTION } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/read_distribution_logs")
include { ASSESS_STRANDEDNESS } from './modules/rseqc.nf'

include { MULTIQC as INFER_EXPERIMENT_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses", MQCLabel:"rseqc-infer_experiment")
include { MULTIQC as GENEBODY_COVERAGE_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses", MQCLabel:"rseqc-genebody_coverage")
include { MULTIQC as INNER_DISTANCE_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses", MQCLabel:"rseqc-inner_distance")
include { MULTIQC as READ_DISTRIBUTION_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses", MQCLabel:"rseqc-read_distribution")


workflow strandedness{
  take:
    bam_array // array: sample-wise tuples (meta, bam_file)
    genome_bed
    samples_ch
 
  main:
     
     bam_array | SORT_INDEX_BAM 
     SORT_INDEX_BAM.out.bam | combine( genome_bed ) |  set { ch_bam_bed }

     INFER_EXPERIMENT( ch_bam_bed )
     GENEBODY_COVERAGE( ch_bam_bed )
     INNER_DISTANCE( ch_bam_bed )
     READ_DISTRIBUTION( ch_bam_bed )
    
     // duplicated in each subworkflow, could use refactoring 
     ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")
 
     INFER_EXPERIMENT_MULTIQC( samples_ch, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config )
     GENEBODY_COVERAGE_MULTIQC( samples_ch, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config )
     INNER_DISTANCE_MULTIQC( samples_ch, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config )
     READ_DISTRIBUTION_MULTIQC( samples_ch, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config )
     
     ch_software_versions = Channel.empty()
     ch_software_versions | mix(INFER_EXPERIMENT.out.version,
                                GENEBODY_COVERAGE.out.version,
                                INNER_DISTANCE.out.version,
                                READ_DISTRIBUTION.out.version)
                          | set{ ch_software_versions }

     INFER_EXPERIMENT.out.log | map { it[1] }
                               | collect
                               | set { ch_infer_expt }
    
     ch_infer_expt | ASSESS_STRANDEDNESS
 
     ch_rseqc_logs = Channel.empty()
     ch_rseqc_logs | mix(INFER_EXPERIMENT.out.log_only,
                         GENEBODY_COVERAGE.out.log_only,
                         INNER_DISTANCE.out.log_only,
                         READ_DISTRIBUTION.out.log_only)
                   | set{ ch_rseqc_logs }
     

  emit:
     strandedness = ASSESS_STRANDEDNESS.out 
     infer_expt = ch_infer_expt
     versions = ch_software_versions 
     rseqc_logs = ch_rseqc_logs
     infer_expt_mqc = INFER_EXPERIMENT_MULTIQC.out.data
}
