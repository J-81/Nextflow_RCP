process THANK_YOU_FOR_COMING_TO_MY_PRESENTATION {
  // conda "${baseDir}/envs/rsem.yml"
  // tag "Sample: ${ meta.id }"
  publishDir "${ params.gldsAccession }"

  input:
    tuple val(meta), path("starOutput/*"), path(RSEM_REF)
  output:

    tuple val(meta), path("questions-${ meta.id }.txt")

  script:
    """
    echo "I'd be happy to take any questions" > questions-${ meta.id }.txt
    """
}
