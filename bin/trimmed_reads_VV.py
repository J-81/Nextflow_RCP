#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from dp_tools.bulkRNASeq.vv_protocols import BulkRNASeq_VVProtocol
from dp_tools.bulkRNASeq.loaders import load_BulkRNASeq_STAGE_00, load_BulkRNASeq_STAGE_01

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser()

    parser.add_argument("--root-path", required=True, help="Root data path")

    parser.add_argument("--accession", required=True, help="Accession number")
    parser.add_argument(
        "--max-flag-code",
        default=80,
        help="Throw an exception if any flag code exceeds this value",
    )

    args = parser.parse_args()
    return args


def main(root_dir: Path, accession: str, max_flag_code: int):
    ds = load_BulkRNASeq_STAGE_01(
        *load_BulkRNASeq_STAGE_00(root_dir, dataSystem_name=accession, stack=True)
    )
    vv_protocol = BulkRNASeq_VVProtocol(
        dataset=ds.dataset, dry_run=False, protocol_name="only trimmed reads"
    )
    vv_protocol.validate_all()
    # output default dataframe
    df = vv_protocol.flags_to_df()
    output_fn = f"VV_log.tsv"
    df.to_csv(output_fn, sep="\t")

    # output verbose dataframe
    df = vv_protocol.flags_to_df(schema="verbose")
    output_fn = f"VV_log_verbose.tsv"
    df.to_csv(output_fn, sep="\t")

    # halt on error
    flagged_messages = '\n'.join([msg for msg in df.loc[df['flag_code'] > max_flag_code]['message']])
    assert (
        df["flag_code"].max() < max_flag_code
    ), f"Maximum flag code exceeded: {max_flag_code}. Printing flag messages that caused this halt: {flagged_messages}"

if __name__ == "__main__":
    args = _parse_args()
    main(
        Path(args.root_path), accession=args.accession, max_flag_code=args.max_flag_code
    )
