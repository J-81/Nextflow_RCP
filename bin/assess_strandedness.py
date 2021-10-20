#! /usr/bin/env python

# a script to assess the dataset-wise strand selection from the sample-wise output from RSeQC infer_experiment
from pathlib import Path
from typing import Tuple
from statistics import mean

def main(root_dir: str):
    results = dict()
    sense_results = list()
    antisense_results = list()
    undetermined_results = list()
    for file in Path(root_dir).glob("*infer_experiment_out"):
        with open(file, "r") as f:
            contents = f.read()
        result = _get_stranded_tuple(contents)
        undetermined_results.append(result[0])
        antisense_results.append(result[1])
        sense_results.append(result[2])

        results[file.name] = result
    
    #print(results)
    # determine average strandedness
    mean_sense  =  mean(sense_results)
    mean_antisense = mean(antisense_results)
    mean_undetermined = mean(undetermined_results)

    if mean_sense > mean_antisense:
        dominant = f"sense:{mean_sense}"
    elif mean_antisense > mean_sense:
        dominant = f"antisense:{mean_antisense}"
    else:
        dominant = f"equal_percents:{mean_sense}"


    with open("result.txt", "w") as f:
        f.write(dominant)
    


def _get_stranded_tuple(text: str) -> Tuple[float,float,float]:
    """ Parses stdout from infer_experiment """
    for line in text.split("\n"):
        if line.startswith("Fraction of reads failed to determine:"):
            undetermined = float(line.split()[-1])
        elif line.startswith('Fraction of reads explained by "1++,1--,2+-,2-+":') or line.startswith('Fraction of reads explained by "++,--":') :
            antisense = float(line.split()[-1])
        elif line.startswith('Fraction of reads explained by "1+-,1-+,2++,2--":') or line.startswith('Fraction of reads explained by "+-,-+":'):
            sense = float(line.split()[-1])
    return (undetermined, sense, antisense)

if __name__ == "__main__":
    main("infer_out")
