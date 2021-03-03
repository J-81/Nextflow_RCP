#! /usr/bin/env python
import zipfile
import tempfile
import os
import sys

from isatools.io import isatab_parser
from isatools.io.isatab_parser import ISATabRecord

def _unzip_ISA(isa_zip_path: str) -> str:
    """ Unzips ISA and places into a tmp contents folder.
    Returns path to temporary directory holding ISA zip file contents.

    :param isa_zip_path: path to isa zip file
    """
    temp_dir = tempfile.mkdtemp()
    with zipfile.ZipFile(isa_zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
    return temp_dir

def _print_metadata(isaobject: ISATabRecord, level: int = 0):
    """ Prints metadata for an ISATabRecord.
    Converts empty strings into None for visibility

    :param isaobject: ISATabRecord as parsed with the isatools library
    :param level: Adds tabs before printing
    """
    for key,metadatum in isaobject.metadata.items():
        # enhance visibility for empty entries
        if (not metadatum and metadatum != 0):
            metadatum = None
            # indent
            print("\t"*level,end='')
            print(f"{key}: {metadatum}")

def parse_isa_dir_from_zip(isa_zip_path: str, pretty_print: bool = False) -> ISATabRecord:
    """ Unzips ISA zip files as found in GLDS metadata folders.

    ISA record metadata signature should match the specs here:
    https://isa-specs.readthedocs.io/en/latest/isamodel.html

    :param isa_zip_path: path to isa zip file
    :param pretty_print: print contents of parsed file, useful for debugging
    """
    INVESTIGATION_INDENT = 1
    STUDY_INDENT = 2
    ASSAY_INDENT = 3

    isa_temp_dir = _unzip_ISA(isa_zip_path)
    investigation = isatab_parser.parse(isatab_ref=isa_temp_dir)

    # only print if requested, useful for debugging parsing and data extraction
    if pretty_print:
        print(f"INVESTIGATION: ISA ZIP: {os.path.basename(isa_zip_path)}")
        print("="*95)
        _print_metadata(investigation, level=INVESTIGATION_INDENT)
        for i, study in enumerate(investigation.studies):
            print("\n")
            print(f"STUDY {i+1} of {len(investigation.studies)}")
            print("="*95)
            _print_metadata(study, level=STUDY_INDENT)
            for design_descriptor in study.design_descriptors:
                [print("\t"*STUDY_INDENT + f"{k}:  {v}") for k,v in design_descriptor.items()]

            # iterature thorugh assays
            for j, assay in enumerate(study.assays):
                print("\n")
                print(f"ASSAY {j+1} of {len(study.assays)} from STUDY {i+1}")
                print("="*95)
                _print_metadata(assay, level=ASSAY_INDENT)

    return investigation

def get_sample_names(isa_zip_path: str,
                     samples_only: bool = False) -> None:
    """ Extracts investigation sample names given a GLDS isa zip file path.
    Returns a dictionary

    :param isa_zip_path: path to isa zip file
    :param samples_only: default, returns dictionary of values, if true, list of sample names only is returned
    """
    samples = dict()
    samples_only_list = list()
    investigation = parse_isa_dir_from_zip(isa_zip_path)
    for study in investigation.studies:
        # study level
        study_key = f"STUDY: {study.metadata['Study Title']}"
        samples[study_key] = dict()
        for assay in study.assays:
            # assay level
            assay_key = f"ASSAY: {assay.metadata['Study Assay Measurement Type']}"
            sample_nodes = [node for node in assay.nodes.values() if node.ntype == "Sample Name"]
            new_samples = [sample_node.name for sample_node in sample_nodes]
            samples[study_key][assay_key] = new_samples
            samples_only_list.extend(new_samples)
    if samples_only:
        return list(set(samples_only_list)) # list,set trick to return only non-redudant set
    return samples

if __name__ == '__main__':
    with open("out_samples.txt", "w") as f:
        samples = get_sample_names(isa_zip_path = sys.argv[1], samples_only = True)
        f.write("\n".join(samples))
