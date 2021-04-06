#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from VV.rsem import RsemCounts
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

    cutoffs = load_cutoffs(None, "DEFAULT")

    cross_checks = dict()
    sample_sheet = RNASeqSampleSheet(sample_sheet = args.runsheet_path)
    cross_checks["SampleSheet"] = sample_sheet
    RsemCounts(dir_mapping = sample_sheet.RSEM_Counts_dir_mapping,
               flagger = flagger,
               has_ERCC = sample_sheet.has_ERCC,
               cutoffs = cutoffs).cross_check
