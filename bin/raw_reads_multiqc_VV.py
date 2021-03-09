#! /usr/bin/env python
""" Validation/Verification for {TODO} in RNASeq Concensus Pipeline
"""
from __future__ import annotations
from collections import namedtuple, defaultdict
import configparser
import argparse
from pathlib import Path
import gzip
import json

from VV.utils import outlier_check, FLAG_LEVELS

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
# TODO: ensure IO signature matches NF process
def _parse_args():
    """ Parse command line args.
    """
    parser = argparse.ArgumentParser(description='Perform Automated V&V on '
                                                 'raw reads for RNASeq Data using '
                                                 'multiQC files')
    parser.add_argument('--config', metavar='c', nargs='+', required=True,
                        help='INI format configuration file')

    parser.add_argument('--samples', required=True,
                        help='Textfile with samples.  Each sample must be in its own line')

    parser.add_argument('--input', required=True,
                        help='Path multiQC data output directory')

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
def validate_verify(samples: list[str], multiQC_data_dir: Path):
    """ Performs VV for raw reads
    """
    # First, data parsing

    # All samples are in one file for multiQC
    file_mapping = defaultdict(lambda: defaultdict(dict))
    file_mapping["all"]["multiQC_json"] = multiQC_data_dir / "multiqc_data.json"

    # TODO: extract data and generate data mapping
    data_mapping = defaultdict(lambda: defaultdict(dict))
    data_mapping = _extract_multiQC_data(file_mapping["all"]["multiQC_json"],
                                         data_mapping,
                                         samples)
    #print(data_mapping)
    print("DATA EXTRACTED FOR EACH SAMPLE")
    [print(key) for key in data_mapping[samples[0]].keys()]
    # create data mapping for 'all' samples entries
    # this is used for outlier detection
    for key in data_mapping[samples[0]].keys():
        data_mapping = _compile_across_samples( data_mapping,
                                                key = key,
                                                samples = samples,
                                                subset_name = "all")
    [print(key, "="*90, value, "\n"*3) for key, value in data_mapping['all'].items()]

    # TODO: perform VV checking
    #   Note: the logging of the VV checks are handled by the write_results function
    vv_mapping = defaultdict(lambda: defaultdict(dict))
    for sample in samples:

        # mapping to save results by sample
        sample_file_mapping = file_mapping[sample]
        sample_data_mapping = data_mapping[sample]
        sample_vv_mapping = vv_mapping[sample]



        sample_vv_mapping["average_read_lengths_match"] = \
            _check_average_read_lengths(sample_data_mapping,
                                        data_mapping)


        if config["GLDS"].getboolean("PairedEnd"):
            sample_vv_mapping["paired_reads_counts_match"] = _check_paired_reads_counts_match(sample_data_mapping)

    vv_mapping["all"]["file_exists"] = (file_mapping["all"]["multiQC_json"].exists(), f"Expected files exist: {file_mapping['all']['multiQC_json'].exists()} "
                                                                        f"File: {file_mapping['all']['multiQC_json']}")

    return file_mapping, data_mapping, vv_mapping

def _check_read_lengths(paired_end, vv_mapping, sample_vv_mapping):
    # TODO: implement with GLDS-104
    if paired_end:
        forward_avg_sequence_length_deviation = \
                outlier_check(value = sample_data_mapping["forward-avg_sequence_length"],
                              against = vv_mapping["all"]["forward-avg_sequence_length"])
        reverse_avg_sequence_length_deviation = \
                outlier_check(value = sample_data_mapping["reverse-avg_sequence_length"],
                              against = vv_mapping["all"]["reverse-avg_sequence_length"])

    else:
        read_avg_sequence_length_deviation = \
                outlier_check(value = sample_data_mapping["read-avg_sequence_length"],
                              against = vv_mapping["all"]["read-avg_sequence_length"])

def _check_average_read_lengths(sample_data_mapping, data_mapping):
    tolerance = config["Raw"].getfloat("SequenceLengthVariationTolerance")
    if  config["GLDS"].getboolean("PairedEnd"):
        forward_deviation = outlier_check(sample_data_mapping["forward-avg_sequence_length"],
                                          data_mapping['all']["forward-avg_sequence_length"])
        reverse_deviation = outlier_check(sample_data_mapping["reverse-avg_sequence_length"],
                                          data_mapping['all']["reverse-avg_sequence_length"])
        if any([forward_deviation > tolerance,
               reverse_deviation > tolerance]):
           result =  (False, f"Average Read length deviation of "
                             f"Forward:{forward_deviation} and Reverse:{reverse_deviation} "
                             f"Greater than allowed ({tolerance})")
        else:
           result =  (True, f"Average Read length deviation of "
                            f"Forward:{forward_deviation} and Reverse:{reverse_deviation} "
                            f"does not exceed {tolerance}")
    else:
        read_deviation = outlier_check(sample_data_mapping["read-avg_sequence_length"],
                                          data_mapping['all']["read-avg_sequence_length"])
        if read_deviation > config["Raw"].getfloat("SequenceLengthVariationTolerance"):
           result =  (False, f"Average Read length deviation of "
                             f"Reads: {read_deviation} "
                             f"Greater than allowed ({tolerance})")
        else:
           result =  (True, f"Average Read length deviation of "
                            f"Reads: {read_deviation} "
                            f"does not exceed {tolerance}")
    return result
def _compile_across_samples(data_mapping, key, samples, subset_name):
    """ Creates an entry for 'all' for value in key.

    Automatically determines if the values are dict like {index: value}
    or a single value per sample.
    """
    aggregate = None
    for sample in samples:
        data = data_mapping[sample][key]
        if type(data) == float:
            # initiate aggregate data type to match
            # if already initiated, this does not affect aggregate
            aggregate = list() if not aggregate else aggregate
            aggregate.append(data)
        elif type(data) == dict:
            # these must be {index:value} dicts of length 1

            # initiate aggregate data type to match
            # if already initiated, this does not affect aggregate
            aggregate = defaultdict(list) if not aggregate else aggregate
            for index, value in data.items():
                aggregate[index].append(value)
        else:
            raise ValueError(f"For {key}, {type(data)} type for data is unexpected.  Aggregation not implemented.")
    # finally add aggregate data to the 'all' entry
    data_mapping[subset_name][key] = aggregate
    return data_mapping

# data extraction functions
# TODO: These should return a value that will be assigned directly to the data mapping
def _extract_multiQC_data(json_file: Path, data_mapping, samples):
    paired_end = config["GLDS"].getboolean("PairedEnd")
    with open(json_file, "r") as f:
        raw_data = json.load(f)

    ###  extract general stats
    for file_data in raw_data["report_general_stats_data"]:
        for file_name, data in file_data.items():
            cur_sample = [sample for sample in samples if sample in file_name][0]
            if paired_end:
                cur_read = "forward" if "_R1" in file_name else "reverse"
            else:
                cur_read = "read"
            for key, value in data.items():
                data_mapping[cur_sample][f"{cur_read}-{key}"] = value

    ### extract plot data

    # extract fastqc_sequence_counts_plot
    for plot_name, data in raw_data["report_plot_data"].items():
        # different kinds of plot data need to be extracted with different methods
        plot_type = data["plot_type"]
        if plot_type == "bar_graph":
            data_mapping = _extract_from_bar_graph(data, plot_name, data_mapping, samples)
        elif plot_type == "xy_line":
            data_mapping = _extract_from_xy_line_graph(data, plot_name, data_mapping, samples)
        else:
            raise ValueError(f"Unknown plot type {plot_type}. Data parsing not implemented for multiQC {plot_type}")

    return data_mapping

def _extract_from_bar_graph(data, plot_name, data_mapping, samples):
    # determine data mapping for samples in multiqc (which are files)
    # and samples in a dataset (which may have a forward and reverse read file)
    paired_end = config["GLDS"].getboolean("PairedEnd")
    mqc_samples_to_samples = dict()
    # this should be a list with one entry
    assert len(data["samples"]) == 1
    for i, mqc_sample in enumerate(data["samples"][0]):
        matching_samples = [sample for sample in samples if sample in mqc_sample]
        # only one sample should map
        assert len(matching_samples) == 1

        sample = matching_samples[0]
        if paired_end and "_R1_" in mqc_sample:
            mqc_samples_to_samples[i] = (sample, "forward")
        elif paired_end and "_R2_" in mqc_sample:
            mqc_samples_to_samples[i] = (sample, "reverse")
        elif not paired_end and "_R1_" in mqc_sample:
            mqc_samples_to_samples[i] = (sample, "read")
        else:
            raise ValueError(
                        f"Unexpected file name format for {mqc_sample} "
                        f"when paired end is {paired_end} and samples are {samples}"
                        )

    # iterate through data from datasets
    # this should be a list with one entry
    assert len(data["datasets"]) == 1
    for sub_data in data["datasets"][0]:
        name = sub_data["name"]
        values = sub_data["data"]
        for i, value in enumerate(values):
            sample, sample_file = mqc_samples_to_samples[i]
            data_mapping[sample][f"{sample_file}-{plot_name}-{name}"] = values[i]
    return data_mapping

def _extract_from_xy_line_graph(data, plot_name, data_mapping, samples):
    # determine data mapping for samples in multiqc (which are files)
    # and samples in a dataset (which may have a forward and reverse read file)
    paired_end = config["GLDS"].getboolean("PairedEnd")

    # Iterate through datasets (each typically representing one sample)
    # Nested list of list, top level list includes percentage and counts

    # plots with multiple value types are detected
    multiple_value_types = True if len(data["datasets"]) > 1 else False

    # plots with bins are detected
    if "categories" in data["config"].keys():
        bins = [str(bin) for bin in data["config"]["categories"]]
    else:
        bins = False


    # dataset represents an entire plot (i.e. all lines)
    # Note: for xy plots with both percent and raw counts, there will be two datasets
    for i, dataset in enumerate(data["datasets"]):
        if multiple_value_types:
            data_label = f"-{data['config']['data_labels'][i]['name']}"
        else:
            data_label = ""

        # dataset entry represents one sample (i.e. one line from the line plot)
        for dataset_entry in dataset:
            file_name = dataset_entry["name"]

            # values is a list of [index,value]
            values = dataset_entry["data"]
            matching_samples = [sample for sample in samples if sample in file_name]
            # only one sample should map
            assert len(matching_samples) == 1

            sample = matching_samples[0]
            if paired_end and "_R1" in file_name:
                sample_file = "forward"
            elif paired_end and "_R2" in file_name:
                sample_file = "reverse"
            elif not paired_end and "_R1" in file_name:
                sample_file = "read"

            # three level nested dict entries for xy graphs
            # {sample: {sample_file-plot_type: {index: value}}}
            data_key = f"{sample_file}-{plot_name}{data_label}"
            data_mapping[sample][data_key] = dict()
            cur_data_mapping_entry = data_mapping[sample][data_key]
            # for plots with bins, add bin string to values iterable
            if bins:
                values = zip(bins, values)
            # for non-categorical bins, each values should be an [index,value]
            for j, value in values:
                cur_data_mapping_entry[j] = value
    return data_mapping

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

def _check_paired_reads_counts_match(sample_data_mapping):
    forward_count = sample_data_mapping["forward-total_sequences"]
    reverse_count = sample_data_mapping["reverse-total_sequences"]
    if forward_count == reverse_count:
        result = (True, f"Paired reads counts match")
    else:
        result = (False, f"Paired reads do not match. Forward count: {forward_count} Reverse count: {reverse_count}")
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
    step_checked = "Raw Reads By MultiQC"
    with open(output, "a+") as f:

        # while less efficient than iterating through samples once, iteration
        #    through samples for each kind of report will make the report more
        #    readable without needing any comment filtering and sorting
        all_vv_mapping = vv_mapping['all']
        key = "file_exists"
        # define report portions that are the same regardless of VV result
        flag_level = None # Overwritten depending on result
        entity = "All Samples"
        details = all_vv_mapping[key][1]
        check_id = f"R_1000"
        check_passed = all_vv_mapping[key][0]

        # Many checks will report both issues and non-issue data
        #   Reminder: entries are tuples with (Pass: bool, Details: str)

        # Issue detected
        if not check_passed:
            flag_level = FLAG_LEVELS[70]
            f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")
        # No Issue detected
        else:
            flag_level = FLAG_LEVELS[20]
            f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")

        for sample in samples:
            sample_vv_mapping = vv_mapping[sample]

            key = "paired_reads_counts_match"
            # define report portions that are the same regardless of VV result
            flag_level = None # Overwritten depending on result
            entity = sample
            details = sample_vv_mapping[key][1]
            check_id = f"R_1001"
            check_passed = sample_vv_mapping[key][0]

            # Many checks will report both issues and non-issue data
            #   Reminder: entries are tuples with (Pass: bool, Details: str)

            # Issue detected
            if not check_passed:
                flag_level = FLAG_LEVELS[70]
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")
            # No Issue detected
            else:
                flag_level = FLAG_LEVELS[20]
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")

        for sample in samples:
            sample_vv_mapping = vv_mapping[sample]

            key = "average_read_lengths_match"
            # define report portions that are the same regardless of VV result
            flag_level = None # Overwritten depending on result
            entity = sample
            details = sample_vv_mapping[key][1]
            check_id = f"R_1002"
            check_passed = sample_vv_mapping[key][0]

            # Many checks will report both issues and non-issue data
            #   Reminder: entries are tuples with (Pass: bool, Details: str)

            # Issue detected
            if not check_passed:
                flag_level = FLAG_LEVELS[60]
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")
            # No Issue detected
            else:
                flag_level = FLAG_LEVELS[20]
                f.write(f"{flag_level}\t{check_id}\t{entity}\t{step_checked}\t{details}\n")

if __name__ == '__main__':
    # Basic comments to add to VV log for each major section
    # Note: each section starts with a double ## line
    # TODO: fill out required logging variables
    VV_name = "Raw Reads (MultiQC)"
    with open(args.output, "a+") as f:
        f.write(f"##{VV_name}\n")
        f.write(f"#{'='*60}\n")
        f.write("#Validation Begin\n")
        for arg_key, arg_value in vars(args).items():
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
    vv_mapping = validate_verify(samples, multiQC_data_dir = Path(args.input))

    # this writes the results to the tsv file
    # the function must be aware of what the vv_mapping will look like
    #   as this will be different for each vv script
    write_results(Path(args.output), vv_mapping)

    # closing line, mainly for visual distinction
    with open(args.output, "a+") as f:
        f.write(f"#{'='*60}\n\n")
