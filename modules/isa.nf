process SAMPLES_FROM_ISA {
  conda "${baseDir}/envs/isatools.yml"

  input:
    path(isaZip)
  output:
    path("out_samples.txt"), emit: samples
  script:
    """
    samples_from_isa.py ${isaZip}
    """
}
