process PARSE_ISA {
    conda "${baseDir}/envs/rnaseq_v1.0_modify.yml"

  input:
    path(isaZip)
  output:
    path("samples.txt"), emit: samples
    path("VV_results.tsv"), emit: vv_results
  script:
    $/
    #! /usr/bin/env python
    from pathlib import Path

    from VV.data import Dataset
    from VV.utils import load_config, load_params
    from VV.flagging import Flagger

    # load parameters and flagger object
    params_file = None if '${ params.vvParamsFile }' == "" else '${ params.vvParamsFile }'
    params_set = None if '${ params.vvParamSet }' == "" else '${ params.vvParamSet }'
    params = load_params(params_file, params_set)
    flagger = Flagger(__file__,
                      halt_level=int('${ params.vvHaltLevel }'),
                      log_to=Path("VV_results.tsv"))

    # Parse ISA into isa object
    isa = Dataset(isa_zip_path = Path('${ isaZip }'),
                  flagger = flagger,
                  entity = '${ params.datasetName }')
    # Perform VV checks.  Uses flagger to log
    isa.validate_verify(vv_for = "RNASeq")

    # Extract samples and save to text file
    samples = isa.get_sample_names(assay_name = "transcription profiling by RNASeq")
    with open("samples.txt", "w") as f:
        for sample in samples:
            f.write(sample + "\n")
    /$
}
