#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from dp_tools.core.post_processing import update_curation_tables
from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
    load_BulkRNASeq_STAGE_02,
    load_BulkRNASeq_STAGE_0201,
    load_BulkRNASeq_STAGE_03,
    load_BulkRNASeq_STAGE_04,
)

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser()

    parser.add_argument("--root-path", required=True, help="Root data path")

    parser.add_argument("--accession", required=True, help="Accession number")

    args = parser.parse_args()
    return args


def main(root_dir: Path, accession: str):
    ds = load_BulkRNASeq_STAGE_04(
        *load_BulkRNASeq_STAGE_03(
            *load_BulkRNASeq_STAGE_0201(
                *load_BulkRNASeq_STAGE_02(
                    *load_BulkRNASeq_STAGE_01(
                        *load_BulkRNASeq_STAGE_00(
                            root_dir, dataSystem_name=accession, stack=True
                        ),
                        stack=True,
                    ),
                    stack=True,
                ),
                stack=True,
            ),
            stack=True,
        )
    )
    update_curation_tables(ds.dataset, config=("bulkRNASeq", "Latest"))


if __name__ == "__main__":
    args = _parse_args()
    main(Path(args.root_path), accession=args.accession)
