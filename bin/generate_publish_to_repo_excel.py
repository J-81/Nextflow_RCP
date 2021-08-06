#! /usr/bin/env python
""" This script generates excel files indicating the files to publish
"""
import sys
from pathlib import Path

import pandas as pd

from displayablepaths import DisplayablePath

def infer_sample(path: Path, sample_list: list) -> str:
    """ Infers the sample based on the file path
    """
    inferred_sample = None
    for sample in sample_list:
        if sample in str(path.path):
            if not inferred_sample:
                inferred_sample = sample
            else:
                raise ValueError(f"Inferred another sample: {sample} but already inferred {inferred_sample}")
    if inferred_sample:
        return inferred_sample
    else:
        #print(f"No sample inferred from {str(path.path)}")
        #print("Returning '?'")
        return 'All samples'


def get_samples(runsheet_path: Path) -> list:
    """ Returns a list of sample names from a given runsheet
    """
    df = pd.read_csv(runsheet_path)
    sample_col = [col for col in df.columns if "sample" in col.lower()][0]
    return list(df[sample_col])


def main(root_path, runsheet_path, template, outputDir: Path = None):
    if template == "RNASeq":
        samples = get_samples(runsheet_path)
        #print(samples)
        paths = DisplayablePath.make_tree(root_path)
        pandas_rows = list()
        for path in paths:
            if path.path.is_file():
                #print(path.path)
                sub_root_dir = path.path.parents[0]
                try:
                    #print(f"{infer_sample(path, samples)}, {path.path.parents[len(path.path.parents)-3].name}, {path.path.name}")
                    
                    parent_dir = path.path.relative_to(root_path).parents
                    print(path.path)
                    print(list(parent_dir))
                    print(list(parent_dir)[0])
                    pandas_rows.append({'Sample Name': infer_sample(path, samples),
                                        'ParentDir':str(list(parent_dir)[0]),
                                        'FileName':path.path.name})
                except IndexError: #files in the root dir GLDS-373/software_versions.txt
                    #print(f"{infer_sample(path, samples)}, {'ROOT'}, {path.path.name}")
                    pandas_rows.append({'Sample Name': infer_sample(path, samples),
                                        'ParentDir':'ROOT',
                                        'FileName':path.path.name})
            else:
                #print(f"Skipping directory: {path}")
                ...

        # create table
        df = pd.DataFrame(pandas_rows).set_index(keys="Sample Name")
        df.to_csv("test.out")

        ########################################################################
        # RNASeq specific conversion to excel tables
        ########################################################################
        ########################################################################
        ###### Raw File Excel File
        ########################################################################
        df_raw_files = df.loc[df["ParentDir"].str.startswith("00-RawData")]\
                         .drop(labels="ParentDir",axis=1)
        Fastq_Files = [file if file.endswith(".fastq.gz") else '' for file in df_raw_files["FileName"] ]
        df_raw_files["FastQ Files"] = Fastq_Files
        df_raw_files.to_csv("test2.out")
        expected_raw_multiqc_filename = "raw_multiqc_report.zip"
        assert expected_raw_multiqc_filename in list(df_raw_files["FileName"]), f"Did not find {expected_raw_multiqc_filename} in 00-RawData"
        df_raw_files["MultiQC Files"] = "raw_multiqc_report.zip"

        # cleanup
        df_raw_files = df_raw_files.drop(labels="FileName", axis=1)
        df_raw_files = df_raw_files.loc[df_raw_files["FastQ Files"] != '']
        df_raw_files = sort_index_and_convert_duplicates(df_raw_files)

        # write to excel file
        output_raw_filename = f"raw_file_names.xlsx"
        print(f"Writing {output_raw_filename}")
        df_raw_files.to_excel(output_raw_filename)
        ########################################################################
        ###### Processed File Excel File
        ########################################################################
        output_processed_name = f"processed_file_names.xlsx"
        with pd.ExcelWriter(output_processed_name) as writer:
            ########################################################################
            ########### 01-TG_Preproc -> Trimmed_Files Tab
            ########################################################################
            df_processed_files = df.loc[df["ParentDir"].str.contains("01-TG_Preproc")]\
                                   .drop(labels="ParentDir",axis=1)
            df_processed_files["Trimming Reports"] = [file if file.endswith(".fastq.gz_trimming_report.txt") else '' for file in df_processed_files["FileName"] ]
            # assert this is still the name for the multiQC file

            # cleanup
            df_processed_files = df_processed_files.drop(labels="FileName", axis=1)
            df_processed_files = df_processed_files.loc[df_processed_files["Trimming Reports"] != '']
            df_processed_files = sort_index_and_convert_duplicates(df_processed_files)

            # write to excel file
            df_processed_files.to_excel(writer, sheet_name="Trimmed_Files")


            ########################################################################
            ########### 02-STAR_Alignment -> Alignment_Files
            ########################################################################
            def ends_with_bam_or_out_tab(filename: str) -> bool:
                if filename.endswith(".bam") or filename.endswith("_SJ.out.tab"):
                    return True
                else:
                    return False

            df_processed_files = df.loc[df["ParentDir"] == "02-STAR_Alignment"]\
                                   .drop(labels="ParentDir",axis=1)
            df_processed_files["Alignment Data"] = [file if ends_with_bam_or_out_tab(file) else '' for file in df_processed_files["FileName"] ]
            # assert this is still the name for the multiQC file

            # cleanup
            df_processed_files = df_processed_files.drop(labels="FileName", axis=1)
            df_processed_files = df_processed_files.loc[df_processed_files["Alignment Data"] != '']
            df_processed_files = sort_index_and_convert_duplicates(df_processed_files)

            # write to excel file
            df_processed_files.to_excel(writer, sheet_name="Alignment_Files")

            ########################################################################
            ########### 03-RSEM_Counts -> Raw_Counts_Files
            ########################################################################
            def allowed_filenames(filename: str) -> bool:
                if filename.endswith(".results") or filename == "Unnormalized_Counts.csv":
                    return True
                else:
                    return False

            df_processed_files = df.loc[df["ParentDir"].isin(["03-RSEM_Counts","04-DESeq2_NormCounts"])]\
                                   .drop(labels="ParentDir",axis=1)
            df_processed_files["Raw Counts Data"] = [file if allowed_filenames(file) else '' for file in df_processed_files["FileName"] ]
            # assert this is still the name for the multiQC file

            # cleanup
            df_processed_files = df_processed_files.drop(labels="FileName", axis=1)
            df_processed_files = df_processed_files.loc[df_processed_files["Raw Counts Data"] != '']
            df_processed_files = sort_index_and_convert_duplicates(df_processed_files)

            # write to excel file
            df_processed_files.to_excel(writer, sheet_name="Raw_Counts_Files")
            ########################################################################
            ########### 04-DESeq2_NormCounts -> Normalized_Counts_Files
            ########################################################################
            def allowed_filenames(filename: str) -> bool:
                if filename == "ERCC_Normalized_Counts.csv" or filename == "Normalized_Counts.csv":
                    return True
                else:
                    return False

            df_processed_files = df.loc[df["ParentDir"].isin(["04-DESeq2_NormCounts"])]\
                                   .drop(labels="ParentDir",axis=1)
            df_processed_files["Normalized Counts Data"] = [file if allowed_filenames(file) else '' for file in df_processed_files["FileName"] ]
            # assert this is still the name for the multiQC file

            # cleanup
            df_processed_files = df_processed_files.drop(labels="FileName", axis=1)
            df_processed_files = df_processed_files.loc[df_processed_files["Normalized Counts Data"] != '']
            df_processed_files = sort_index_and_convert_duplicates(df_processed_files)

            # write to excel file
            df_processed_files.to_excel(writer, sheet_name="Normalized_Counts_Files")
            ########################################################################
            ########### 05-DESeq2_DGE -> DGE_Files
            ########################################################################
            def allowed_filenames(filename: str) -> bool:
                if filename in ["contrasts.csv","differential_expression.csv",
                "visualization_output_table.csv","visualization_PCA_table.csv",
                "ERCCnorm_contrasts.csv","ERCCnorm_differential_expression.csv",
                "visualization_output_table_ERCCnorm.csv","visualization_PCA_table_ERCCnorm.csv"]:
                    return True
                else:
                    return False

            df_processed_files = df.loc[df["ParentDir"].isin(["05-DESeq2_DGE"])]\
                                   .drop(labels="ParentDir",axis=1)
            df_processed_files["DGE Data"] = [file if allowed_filenames(file) else '' for file in df_processed_files["FileName"] ]
            # assert this is still the name for the multiQC file

            # cleanup
            df_processed_files = df_processed_files.drop(labels="FileName", axis=1)
            df_processed_files = df_processed_files.loc[df_processed_files["DGE Data"] != '']
            df_processed_files = sort_index_and_convert_duplicates(df_processed_files)

            # write to excel file
            print(f"Writing {output_processed_name}")
            df_processed_files.to_excel(writer, sheet_name="DGE_Files")



def sort_index_and_convert_duplicates(df, duplicateTo=''):
    already_used_index = list()
    df = df.sort_index()
    copy_index = df.index.copy()
    new_index = list()
    for i, index in enumerate(copy_index):
        if index in already_used_index:
            new_index.append('')
        else:
            new_index.append(index)
        already_used_index.append(index)
    df.index = pd.Index(new_index, name=copy_index.name)
    return df

if __name__ == '__main__':
    main(root_path = Path(sys.argv[1]), runsheet_path = Path(sys.argv[2]), template = "RNASeq")
