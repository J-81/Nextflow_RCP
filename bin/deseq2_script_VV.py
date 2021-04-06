#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path
import tempfile

from VV.rsem import RsemCounts
from VV.deseq2 import Deseq2ScriptOutput
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
#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from VV import raw_reads
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

    null_flagger = Flagger(script = __file__,
                           log_to = Path(tempfile.TemporaryFile()),
                           halt_level = int(args.halt_severity))

    cutoffs = load_cutoffs(None, "DEFAULT")

    cross_checks = dict()
    sample_sheet = RNASeqSampleSheet(sample_sheet = args.runsheet_path)
    cross_checks["SampleSheet"] = sample_sheet

    rsem_cross_check =   RsemCounts(dir_mapping = sample_sheet.RSEM_Counts_dir_mapping,
                                    flagger = null_flagger, # we don't want to reflag
                                    has_ERCC = sample_sheet.has_ERCC,
                                    cutoffs = cutoffs).cross_check
    cross_checks["RSEM"] = rsem_cross_check

    Deseq2ScriptOutput(samples = sample_sheet.samples,
                       counts_dir_path = sample_sheet.DESeq2_NormCount,
                       dge_dir_path = sample_sheet.DESeq2_DGE,
                       flagger = flagger,
                       cutoffs = cutoffs,
                       has_ERCC = sample_sheet.has_ERCC,
                       cross_checks = cross_checks)
    ###########################################################################
    # Generate derivative log files
    ###########################################################################
    print(f"{'='*40}")
    for log_type in ["only-issues", "by-sample", "by-step","all-by-entity"]:
        flagger.generate_derivative_log(log_type = log_type,
                                        samples = sample_sheet.samples)
