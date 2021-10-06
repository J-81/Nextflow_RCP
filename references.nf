// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2

include { DOWNLOAD_ERCC } from './modules/download.nf'
include { CONCAT_ERCC;
          SUBSAMPLE_GENOME;
          TO_PRED;
          TO_BED } from './modules/genome.nf'

include { DOWNLOAD_GENOME_ANNOTATIONS as DOWNLOAD_TOPLEVEL_REF } from './modules/download.nf' addParams(ref_target: "toplevel", _has_fallback: false)

include { DOWNLOAD_GENOME_ANNOTATIONS as DOWNLOAD_PRIMARY_ASSEMBLY_REF } from './modules/download.nf' addParams(ref_target: "primary_assembly", _has_fallback: true)

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow references{
  take:
    organism_sci
    has_ercc
  main:



      if ( params.ref_order == 'primary_assemblyELSEtoplevel' ) {
        annotations = Channel.empty()
        DOWNLOAD_PRIMARY_ASSEMBLY_REF( organism_sci ) | map { it -> [2, it] } // add a priority value
                                                      | set{pa}
        DOWNLOAD_TOPLEVEL_REF( organism_sci ) | map {it -> [1, it] } // add a priority value
                                              | set{tl}
        pa.ifEmpty([0,null]) | set {pa}
        
        annotations | mix(pa,tl)
                    | view
                    | max { it[0] }
                    | map { it[1] }
                    | set { genome_annotations_pre_subsample }
        //
        /* 
        annotations | max { it[0] }  // take the highest priority that lands in the completed channel
                    | set { genome_annotations_pre_subsample }
        */
      } else if (params.ref_order == 'toplevel' ) {
      	DOWNLOAD_TOPLEVEL_REF( organism_sci ) | set { genome_annotations_pre_subsample }
      }

      // SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
      if ( params.genomeSubsample ) {
        SUBSAMPLE_GENOME( genome_annotations_pre_subsample, organism_sci )
        SUBSAMPLE_GENOME.out.build | first | set { genome_annotations_pre_ercc }
      } else {
        genome_annotations_pre_subsample | first | set { genome_annotations_pre_ercc }
      }

      // ERCC STEP : ADD ERCC Fasta and GTF to genome files
      CONCAT_ERCC( genome_annotations_pre_ercc, DOWNLOAD_ERCC(), organism_sci, has_ercc )
      .ifEmpty { genome_annotations_pre_ercc.value }  | set { genome_annotations }


      TO_PRED( genome_annotations | map { it[1] } )
      TO_BED( TO_PRED.out, organism_sci )

  emit:
      genome_annotations = genome_annotations
      genome_bed = TO_BED.out
}
