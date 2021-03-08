#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
from __future__ import annotations
import datetime
from collections import namedtuple, defaultdict
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

from VV.utils import outlier_check, FLAG_LEVELS

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

    parser.add_argument('--samples', required=True,
                        help='Samples to Process')

    parser.add_argument('--input', required=True,
                        help='Paths to raw read files for')

    parser.add_argument('--output', metavar='o', required=True,
                        help='File to write VV results to')


    args = parser.parse_args()
    print(args)
    return args


args = _parse_args()
config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
config.read(args.config)

def validate_verify(samples: list[str], raw_reads_dir: Path):
    """ Performs VV for raw reads
    """
    # First, data parsing

    # generate sample to read mapping
    file_mapping = dict()
    for sample in samples:

        # set up each sample entry as a dictionary
        file_mapping[sample] = dict()

        if config["GLDS"].getboolean("PairedEnd"):
            # file patterns for paired end studies
            # note: this may be replaced in the future using expected filenames specified in the ISA
            file_mapping[sample]["forward_read"] = raw_reads_dir / f"{sample}_R1_raw.fastq.gz"
            file_mapping[sample]["reverse_read"] = raw_reads_dir / f"{sample}_R2_raw.fastq.gz"

        else:
            # file pattern expected for single end studies
            # note: this may be replaced in the future using expected filenames specified in the ISA
            file_mapping[sample]["read"] = raw_reads_dir / f"{sample}_R1_raw.fastq.gz"

    vv_mapping = defaultdict(dict)

    file_sizes = list()
    for sample in samples:
        if config["GLDS"].getboolean("PairedEnd"):
            file_size = file_mapping[sample]["forward_read"].stat().st_size + file_mapping[sample]["reverse_read"].stat().st_size
        else:
            file_size = file_mapping[sample]["read"].stat().st_size + file_mapping[sample]["reverse_read"].stat().st_size
        file_sizes.append(file_size)
        vv_mapping[sample]["file_size"] = file_size
    vv_mapping["all"]["file_sizes"] = file_sizes

    # sampleWise data checks
    for sample in samples:

        # mapping to save results by sample
        sample_file_mapping = file_mapping[sample]
        sample_vv_mapping = vv_mapping[sample]

        sample_vv_mapping["file_exists"] = _check_file_existence(sample_file_mapping)
        sample_vv_mapping["header_check"] = _check_headers(sample_file_mapping.values(), config["Options"].getint("MaxFastQLinesToCheck"))
        sample_vv_mapping["file_size_deviation"] = outlier_check(value = sample_vv_mapping["file_size"],
                                                                 against = vv_mapping["all"]["file_sizes"]
                                                                )

    return file_mapping, vv_mapping

def _check_file_existence(sample_file_mapping):
    missing_files = list()
    for file_alias, file in sample_file_mapping.items():
        if not file.exists():
            missing_files.append(file)
    if len(missing_files) != 0:
        result = (False, f"{missing_files} missing")
    else:
        result = (True, f"All expected files present: {sample_file_mapping.items()}")
    return result

def _check_headers(files, count_lines_to_check: int) -> int:
    """ Checks fastq lines for expected header content

    Note: Example of header from GLDS-194

    |  ``@J00113:376:HMJMYBBXX:3:1101:26666:1244 1:N:0:NCGCTCGA\n``

    This also assumes the fastq file does NOT split sequence or quality lines
    for any read

    :param file: compressed fastq file to check
    :param count_lines_to_check: number of lines to check. Special value: -1 means no limit, check all lines.
    """
    if count_lines_to_check == -1:
        count_lines_to_check = float("inf")

    # TODO: add expected length check
    expected_length = None

    lines_with_issues = list()

    passes = True
    message = ""
    for file in files:
        with gzip.open(file, "rb") as f:
            for i, line in enumerate(f):
                # checks if lines counted equals the limit input
                if i+1 == count_lines_to_check:
                    print(f"Reached {count_lines_to_check} lines, ending line check")
                    break

                line = line.decode()
                # every fourth line should be an identifier
                expected_identifier_line = (i % 4 == 0)
                # check if line is actually an identifier line
                if (expected_identifier_line and line[0] != "@"):
                    lines_with_issues.append(i+1)
                    print(f"FAIL: {checkname}: "
                          f"Line {i+1} of {file} was not an identifier line as expected "
                          f"LINE {i+1}: {line}")
                # update every 20,000,000 reads
                if i % 20000000 == 0:
                    print(f"Checked {i} lines for {file}")
        if len(lines_with_issues) != 0:
            passes = False
            message += f"for {file}, first ten lines with header issues: {lines_with_issues[0:10]} of {len(lines_with_issues)} header lines with issues: "
        else:
            message += f"for {file}, No issues with headers checked up to line {count_lines_to_check}: "
    return (passes, message)

def write_results(output: Path, all_results):
    """ Write VV results to output file.
    Also formats for consistency across VV steps
    """
    with open(output, "a+") as f:
        for sample in samples:
            results = all_results[sample]

            # header check logging
            if results["header_check"][0] == False:
                flag_level = FLAG_LEVELS[50]
                entity = sample
                file_checked = "Raw Reads"
                details = results['header_check'][1]
                check_id = f"R_0001"
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{file_checked}\t{details}\n")
            else:
                flag_level = FLAG_LEVELS[20]
                entity = sample
                file_checked = "Raw Reads"
                details = results['header_check'][1]
                check_id = f"R_0001"
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{file_checked}\t{details}\n")

        for sample in samples:
            results = all_results[sample]

            # file existence logging
            if results["file_exists"][0] == False:
                flag_level = FLAG_LEVELS[70]
                entity = sample
                file_checked = "Raw Reads"
                details = results['file_exists'][1]
                check_id = f"X_0001"
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{file_checked}\t{details}\n")
            else:
                flag_level = FLAG_LEVELS[20]
                entity = sample
                file_checked = "Raw Reads"
                details = results['file_exists'][1]
                check_id = f"X_0001"
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{file_checked}\t{details}\n")

        # file size deviation
        for sample in samples:
            deviation = all_results[sample]["file_size_deviation"]
            flag_level = None
            entity = sample
            file_checked = "Raw Reads"
            details = None
            check_id = f"R_0002"

            if deviation > config["Raw"].getfloat("FileSizeVariationToleranceRed"):
                flag_level = FLAG_LEVELS[50]
                details = f"File size(s) differ from median by {deviation:.2f}"
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{file_checked}\t{details}\n")

            elif deviation > config["Raw"].getfloat("FileSizeVariationToleranceYellow"):
                flag_level = FLAG_LEVELS[50]
                details = f"File size(s) differ from median by {deviation:.2f}"
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{file_checked}\t{details}\n")




if __name__ == '__main__':
    with open(args.output, "a+") as f:
        f.write("##Raw_Reads\n")
        f.write(f"#{'='*60}\n")
        f.write("#Validation Begin\n")
        f.write(f"#Config: {type(args.config)} {args.config}\n")
        f.write(f"#Samples: {type(args.samples)} {args.samples}\n")
        f.write(f"#Input: {type(args.input)} {args.input}\n")
        f.write(f"#Output: {type(args.output)} {args.output}\n")
        f.write("#Validation Results Below\n")


    with open(args.samples, "r") as f:
        samples = [sample.strip() for sample in f.readlines()]

    _, vv_mapping = validate_verify(samples, Path(args.input))

    write_results(Path(args.output), vv_mapping)
    with open(args.output, "a+") as f:
        f.write(f"#{'='*60}\n\n")
