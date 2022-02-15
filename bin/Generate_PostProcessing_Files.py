#! /usr/bin/env python
import sys
import logging

logging.basicConfig(level=logging.DEBUG)

from gl4gl.PathAnnotate.file_search import get_annotated_paths_df

from gl4gl.PostProcessing import generate_files_reports

from gl4gl import PathAnnotate
from pathlib import Path

config_fs_f = [
    c for c in PathAnnotate.get_configs() if c.name == "Bulk_Search_Patterns.yaml"
][0]

template = [
    t for t in PathAnnotate.get_templates(config_fs_f) if t == "Bulk_RNASeq:PairedEnd"
][0]


def main():
    processed_path, raw_path, unpublished_path = generate_files_reports(
        root_path=Path(sys.argv[1]),
        runsheet_path=Path(sys.argv[2]),
        template=template,
        config_fs_f=config_fs_f,
    )

    assert processed_path.name == "processed_file_names.xlsx"
    assert raw_path.name == "raw_file_names.xlsx"
    assert unpublished_path.name == "unpublished_file_names.xlsx"


if __name__ == "__main__":
    main()
