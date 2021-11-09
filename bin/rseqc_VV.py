#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from VV.rseqc import Rseqc
from VV.utils import load_cutoffs
from VV.flagging import Flagger
from VV.rnaseq_samplesheet import RNASeqSampleSheet

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
def _parse_args():
    """ Parse command line args.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--runsheet-path', required=True,
                        help='run sheet path')

    parser.add_argument('--output', metavar='o', required=True,
                        help='File to write VV results to')

    parser.add_argument('--halt-severity', metavar='n', required=True,
                        help='Flag Level to raise an error and halt processing')


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = _parse_args()

    flagger = Flagger(script = __file__,
                      log_to = Path(args.output),
                      halt_level = int(args.halt_severity))

    cutoffs = load_cutoffs(None, "DEFAULT_RNASEQ")

    cross_checks = dict()
    sample_sheet = RNASeqSampleSheet(sample_sheet = args.runsheet_path)
    Rseqc(multiqc_json = sample_sheet.rseqc_multiqc,
          samples = sample_sheet.samples,
          flagger = flagger,
          cutoffs = cutoffs)