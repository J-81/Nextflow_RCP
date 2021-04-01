#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import datetime
from collections import namedtuple
from typing import Tuple
import hashlib
import statistics
import glob
import os
import sys
import configparser
import argparse
from pathlib import Path
import gzip

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
def _parse_args():
    """ Parse command line args.
    """
    parser = argparse.ArgumentParser(description='Perform Automated V&V on '
                                                 'raw and processed RNASeq Data.')
    parser.add_argument('--config', metavar='c', nargs='+', required=True,
                        help='INI format configuration file')

    parser.add_argument('--samples', metavar='s', nargs='+', required=True,
                        help='Samples to Process')

    parser.add_argument('--geneCounts', metavar='g', nargs='+', required=True,
                        help='Paths to .gene.results file created by RSEM')

    parser.add_argument('--transcriptCounts', metavar='t', nargs='+', required=True,
                        help='Paths to .isoform.results file created by RSEM')

    parser.add_argument('--stats', metavar='l', nargs='+', required=True,
                        help='Stat directory from RSEM')


    parser.add_argument('--output', metavar='o', required=True,
                        help='File to write VV results to')


    args = parser.parse_args()
    print(args)
    return args


args = _parse_args()
config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
config.read(args.config)

if __name__ == '__main__':
    with open(args.output, "a+") as f:
        f.write("RSEM\n")
        f.write(f"{'='*60}\n")
        f.write("Validation Begin\n")
        f.write(f"Config: {type(args.config)} {args.config}\n")
        f.write(f"Samples: {type(args.samples)} {args.samples}\n")
        f.write(f"Input_geneCounts: {type(args.geneCounts)} {args.geneCounts}\n")
        f.write(f"Input_transcriptCounts: {type(args.transcriptCounts)} {args.transcriptCounts}\n")
        f.write(f"Input_stats: {type(args.stats)} {args.stats}\n")
        f.write(f"Output: {type(args.output)} {args.output}\n")
        f.write("Validation End\n")
        f.write(f"{'='*60}\n\n")
