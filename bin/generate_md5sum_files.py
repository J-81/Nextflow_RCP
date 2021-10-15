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
    fname = row.filename
    if any(
            (
              (fname.endswith("_raw.fastq.gz")),
              (fname == "raw_multiqc_report.zip"),
              # processed
                # trimmed files
                (fname.endswith("_trimming_report.txt")),
                (fname.endswith("_trimmed.fastq.gz")),
                (fname == "trimmed_multiqc_report.zip"),
                # alignment files
                (fname.endswith("_Aligned.sortedByCoord.out.bam")),
                (fname.endswith("_Aligned.toTranscriptome.out.bam")),
                (fname.endswith("_SJ.out.tab")),

                (fname.endswith("_Log.final.out")),
                (fname == "align_multiqc_report.zip"),

                # raw counts files
                (fname.endswith(".genes.results")),
                (fname.endswith(".isoforms.results")),
                (fname == "Unnormalized_Counts.csv"),

                # Normalized Counts
                (fname == "Normalized_Counts.csv"),
                (fname == "ERCC_Normalized_Counts.csv"),

                # DGE Files
                (fname == "contrasts.csv"),
                (fname == "differential_expression.csv"),
                (fname == "visualization_output_table.csv"),
                (fname == "visualization_PCA_table.csv"),
                (fname == "ERCCnorm_contrasts.csv"),
                (fname == "ERCCnorm_differential_expression.csv"),
                (fname == "visualization_output_table_ERCCnorm.csv"),
                (fname == "visualization_PCA_table_ERCCnorm.csv"),
              )
            ):
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
    md_table.to_csv(output_table_filename.replace(".md",".tsv"), index=False, sep="\t")
    # repeat writing after filtering files the won't be published to the repo
    md_table = md_table.loc[md_table.apply(filter_to_publish_files, axis=1)]
    #md_table =
    output_table_filename = f"{outputDir}/{root_path.name}_md5sum_table_publish_to_repo.md" if outputDir else f"{root_path.name}_md5sum_table_publish_to_repo.md"
    print(f"Writing {output_table_filename}")
    md_table.to_markdown(buf=output_table_filename, index=False)
    md_table.to_csv(output_table_filename.replace(".md",".tsv"), index=False, sep="\t")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
