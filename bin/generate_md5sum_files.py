#! /usr/bin/env python
""" This scripts generates a file tree with md5sums for each file computed

Based on original code sourced from:
https://stackoverflow.com/questions/9727673/list-directory-tree-structure-in-python
- Original code Author: https://stackoverflow.com/users/2479038/abstrus

"""
import sys
from pathlib import Path

import pandas as pd

from displayablepaths import DisplayablePath

def filter_to_publish_files(row):
    fn = row.filename
    if any(
            (
                fn.endswith(".fastq.gz"),
                fn.endswith(".fastq.gz_trimming_report.txt"),
                fn.endswith(".out.bam"),
                fn.endswith("_SJ.out.tab"),
                fn.endswith(".genes.results"),
                fn.endswith(".isoforms.results"),
                fn == "Unnormalized_Counts.csv",
                fn == "ERCC_Normalized_Counts.csv",
                fn == "Normalized_Counts.csv",
                fn == "contrasts.csv",
                fn == "differential_expression.csv",
                fn == "ERCCnorm_contrasts.csv",
                fn == "ERCCnorm_differential_expression.csv",
                fn == "visualization_output_table.csv",
                fn == "visualization_PCA_table.csv",
                fn == "visualization_output_table_ERCCnorm.csv",
                fn == "visualization_PCA_table_ERCCnorm.csv",
            )
        ) :
        return True
    return False

def main(root_path: Path, outputDir: Path = None):
    paths = DisplayablePath.make_tree(root_path)
    output_tree_filename = f"{outputDir}/{root_path.name}_md5sum_tree.txt" if outputDir else f"{root_path.name}_md5sum_tree.txt"
    print(f"Writing {output_tree_filename}")
    with Path(output_tree_filename).open("w") as tree_out:
        pandas_rows = list()
        for path in paths:
            tree_out.write(path.displayable()+'\n')
            if path.path.is_file():
                pandas_rows.append(path.pandas_row)
    md_table = pd.DataFrame(pandas_rows)
    output_table_filename = f"{outputDir}/{root_path.name}_md5sum_table.md" if outputDir else f"{root_path.name}_md5sum_table.md"
    print(f"Writing {output_table_filename}")
    md_table.to_markdown(buf=output_table_filename, index=False)
    # repeat writing after filtering files the won't be published to the repo
    md_table = md_table.loc[md_table.apply(filter_to_publish_files, axis=1)]
    #md_table =
    output_table_filename = f"{outputDir}/{root_path.name}_md5sum_table_publish_to_repo.md" if outputDir else f"{root_path.name}_md5sum_table_publish_to_repo.md"
    print(f"Writing {output_table_filename}")
    md_table.to_markdown(buf=output_table_filename, index=False)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
