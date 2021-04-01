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

    parser.add_argument('--normDir', required=True,
                        help='Path to directory with normalized counts output')

    parser.add_argument('--dgeDir', required=True,
                        help='Path to directory with dge analysis output')

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
        f.write("Deseq2 Script\n")
        f.write(f"{'='*60}\n")
        f.write("Validation Begin\n")
        f.write(f"Config: {type(args.config)} {args.config}\n")
        f.write(f"Samples: {type(args.samples)} {args.samples}\n")
        f.write(f"normDir: {type(args.normDir)} {args.normDir}\n")
        f.write(f"dgeDir: {type(args.dgeDir)} {args.dgeDir}\n")
        f.write(f"Output: {type(args.output)} {args.output}\n")
        f.write("Validation End\n")
        f.write(f"{'='*60}\n\n")
