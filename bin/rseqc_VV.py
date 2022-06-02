#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from dp_tools.bulkRNASeq.vv_protocols import validate_bulkRNASeq
from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
    load_BulkRNASeq_STAGE_02,
    load_BulkRNASeq_STAGE_0201,
)

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
    ds = load_BulkRNASeq_STAGE_0201(
        *load_BulkRNASeq_STAGE_02(
            *load_BulkRNASeq_STAGE_01(
                *load_BulkRNASeq_STAGE_00(
                    root_dir, dataSystem_name=accession, stack=True
                ),
                stack=True,
            ),
            stack=True,
        )
    )
    vp = validate_bulkRNASeq(
        ds.dataset,
        report_args={"include_skipped": True},
        protocol_args={
            "skip_components": [
                "Metadata",
                "Raw Reads",
                "Raw Reads By Sample",
                "Trim Reads",
                "Trimmed Reads By Sample",
                "STAR Alignments",
                "STAR Alignments By Sample",
                # "RSeQC By Sample",
                # "RSeQC",
                "RSEM Counts",
                "Unnormalized Gene Counts",
                "DGE Metadata",
                "DGE Metadata ERCC",
                "DGE Output",
                "DGE Output ERCC",
            ]
        },
        defer_run=True,
    )

    print(f"{'QUEUED CHECK COMPONENT TREE':*^60}")
    print(vp.queued_checks(include_individual_checks=False))
    print(f"{'QUEUED CHECK COMPONENT TREE WITH SKIPPED COMPONENTS':*^60}")
    print(
        vp.queued_checks(
            include_individual_checks=False, include_skipped_components=True
        )
    )
    print(f"{'QUEUED CHECK COMPONENT TREE WITH INVIDUAL CHECKS':*^60}")
    print(vp.queued_checks(include_individual_checks=True))

    vp.run()
    report = vp.report(include_skipped=False)

    # output default dataframe
    df = report["flag_table"]
    output_fn = f"VV_log.tsv"
    df.to_csv(output_fn, sep="\t")

    # halt on error
    flagged_messages = "\n".join(
        [msg for msg in df.loc[df["code_level"] >= max_flag_code]["message"]]
    )
    assert (
        df["code_level"].max() < max_flag_code
    ), f"Maximum flag code exceeded: {max_flag_code}. Printing flag messages that caused this halt: {flagged_messages}"


if __name__ == "__main__":
    args = _parse_args()
    main(
        Path(args.root_path),
        accession=args.accession,
        max_flag_code=int(args.max_flag_code),
    )
