#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from dp_tools.bulkRNASeq.loaders import (
    load_BulkRNASeq_STAGE_00,
    load_BulkRNASeq_STAGE_01,
    load_BulkRNASeq_STAGE_02,
    load_BulkRNASeq_STAGE_0201,
    load_BulkRNASeq_STAGE_03,
    load_BulkRNASeq_STAGE_04,
)
from dp_tools.core.post_processing import (
    ALLOWED_MISSING_KEYS_FOR_NON_ERCC,
    ALLOWED_MISSING_KEYS_FOR_PAIRED_END,
    ALLOWED_MISSING_KEYS_FOR_SINGLE_END,
    generate_md5sum_table,
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

    # determine what data assets can be missing
    missing_keys: set[str] = set()
    if ds.dataset.metadata.paired_end:
        missing_keys = missing_keys.union(ALLOWED_MISSING_KEYS_FOR_PAIRED_END)
    else:
        missing_keys = missing_keys.union(ALLOWED_MISSING_KEYS_FOR_SINGLE_END)
    if not ds.dataset.metadata.has_ercc:
        missing_keys = missing_keys.union(ALLOWED_MISSING_KEYS_FOR_NON_ERCC)

    df = generate_md5sum_table(
        ds.dataset,
        config=("bulkRNASeq", "Latest"),
        allowed_unused_keys=missing_keys,
        include_tags=True,
    )

    unique_tags = set(df["tags"].sum())
    for tag in unique_tags:
        df_subset = df.loc[df["tags"].apply(lambda l: tag in l)].drop(
            "tags", axis="columns"
        )
        df_subset.to_csv(f"{accession}_{tag}_md5sum.tsv", sep="\t", index=False)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    args = _parse_args()
    main(Path(args.root_path), accession=args.accession)
