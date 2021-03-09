#! /usr/bin/env python
""" Validation/Verification for {TODO} in RNASeq Concensus Pipeline
"""
from __future__ import annotations
from collections import namedtuple, defaultdict
import configparser
import argparse
from pathlib import Path
import gzip

from VV.utils import outlier_check, FLAG_LEVELS

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
# TODO: ensure IO signature matches NF process
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
    return args


args = _parse_args()
config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
config.read(args.config)

# TODO: Using Inputs, perform data extraction.
# Three kinds of mappings are utilized for most VV
# file_mapping: {sample: {file_alias: path}}
# data_mapping: {sample: {data_name: data}, all: {data_name: data}}
#     note: all is a special key for storing data about all samples including aggregate data
# vv_mapping: {sample: {check_name: results}, other_entity: {check_name: results}}
#     note: other_entity MUST be tagged with "NS_{entity}" otherwise it will be assumed to be a sample name
def validate_verify(samples: list[str], raw_reads_dir: Path):
    """ Performs VV for raw reads
    """
    # First, data parsing

    # TODO: generate file mapping object
    file_mapping = defaultdict(dict)
    for sample in samples:
        if config["GLDS"].getboolean("PairedEnd"):
            # file patterns for paired end studies
            # note: this may be replaced in the future using expected filenames specified in the ISA
            file_mapping[sample]["forward_read"] = raw_reads_dir / f"{sample}_R1_raw.fastq.gz"
            file_mapping[sample]["reverse_read"] = raw_reads_dir / f"{sample}_R2_raw.fastq.gz"

        else:
            # file pattern expected for single end studies
            # note: this may be replaced in the future using expected filenames specified in the ISA
            file_mapping[sample]["read"] = raw_reads_dir / f"{sample}_R1_raw.fastq.gz"

    # TODO: extract data and generate data mapping
    data_mapping = defaultdict(dict)


    # TODO: perform VV checking
    #   Note: the logging of the VV checks are handled by the write_results function
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

# data extraction functions
# TODO: These should return a value that will be assigned directly to the data mapping
def _extract_data(sample_file_mapping):
    value = None
    return value

# vv check functions
# TODO: These should return a (True/False, String) that will be interpretted by
#   the write results function
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

# writes vv_mapping
# TODO: Each line MUST conform to the example format to ensure the TSV file
#   can be parsed and analysed downstream in a consistent manner
#   Format: {flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n
def write_results(output: Path, vv_mapping):
    """ Write VV results to output file.
    Also formats for consistency across VV steps
    """
    # The step checked will be the same for all reports in a given vv script
    step_checked = "Raw Reads"
    with open(output, "a+") as f:

        # while less efficient than iterating through samples once, iteration
        #    through samples for each kind of report will make the report more
        #    readable without needing any comment filtering and sorting
        for sample in samples:
            results = vv_mapping[sample]

            # define report portions that are the same regardless of VV result
            vv_key = 'TODO'
            flag_level = None # Overwritten depending on result
            entity = sample
            details = results[vv_key][1]
            check_id = f"R_0001"
            check_passed = results[vv_key][0]

            # Many checks will report both issues and non-issue data
            #   Reminder: entries are tuples with (Pass: bool, Details: str)

            # Issue detected
            if not check_passed:
                flag_level = FLAG_LEVELS[50]
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")
            # No Issue detected
            else:
                flag_level = FLAG_LEVELS[20]
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")

if __name__ == '__main__':
    # Basic comments to add to VV log for each major section
    # Note: each section starts with a double ## line
    # TODO: fill out required logging variables
    VV_name = "ADD_HERE"
    with open(args.output, "a+") as f:
        f.write(f"##{VV_name}\n")
        f.write(f"#{'='*60}\n")
        f.write("#Validation Begin\n")
        for arg_key, arg_value in args.items():
            f.write(f"#{arg_key}: {type(arg_value)} {arg_value}\n")
        f.write("#Validation Results Below\n")

    # read in sample names
    with open(args.samples, "r") as f:
        samples = [sample.strip() for sample in f.readlines()]

    # run validate_verify function
    # this performs data extraction and vv checks
    # it does not handle reporting
    # this is to ensure the analysis and reporting of results are decoupled
    # TODO: ensure input signature matches function
    vv_mapping = validate_verify(samples, Path(args.input))

    # this writes the results to the tsv file
    # the function must be aware of what the vv_mapping will look like
    #   as this will be different for each vv script
    write_results(Path(args.output), vv_mapping)

    # closing line, mainly for visual distinction
    with open(args.output, "a+") as f:
        f.write(f"#{'='*60}\n\n")
