#! /usr/bin/env python
import sys
import hashlib
import logging
from pathlib import Path

logging.basicConfig(level=logging.DEBUG)

import pandas as pd

from gl4gl import PathAnnotate


def get_filesize(path: Path):
    """ Used as a test function since md5sum is expensive"""
    logging.debug(f"Getting file size in bytes for: {str(path)}")
    return path.stat().st_size


def get_md5sum(path: Path):
    logging.debug(f"Getting md5sum for: {str(path)}")
    with open(path, "rb") as f:
        md5sum = hashlib.md5(f.read()).hexdigest()
    return md5sum


def main(root_dir: Path, runsheet: Path, compute_attr: str):
    def filter_function(df) -> pd.DataFrame:
        return df.loc[~(df["isDir"] | df["Excel File"] == "unpublished")]

    config_fs_f = [
        c for c in PathAnnotate.get_configs() if c.name == "Bulk_Search_Patterns.yaml"
    ][0]

    template = [
        t
        for t in PathAnnotate.get_templates(config_fs_f)
        if t == "Bulk_RNASeq:PairedEnd"
    ][0]

    compute_func = {"filesize": get_filesize, "md5sum": get_md5sum}[compute_attr]

    compute_these = (
        filter_function,
        [((compute_func, {}), (compute_attr),),],
    )

    df = PathAnnotate.get_annotated_paths_df(
        root_dir=root_dir,
        runsheet=runsheet,
        template=template,
        config_fs_f=config_fs_f,
        drop_annotations=True,
        compute_these=compute_these,
    )

    df["filename"] = df["pathObj"].apply(lambda p: p.name)
    df[["filename", compute_attr]].to_csv("md5sum.tsv", sep="\t", index=None)


if __name__ == "__main__":
    main(
        root_dir=Path(sys.argv[1]),
        runsheet=Path(sys.argv[2]),
        compute_attr=sys.argv[3],
    )
