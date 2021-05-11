# NASA GeneLab Pipeline Reference: GL-DPPD-7101-C

This is a nextflow implementation of the GeneLab RNASeq Concensus Pipeline.

This pipeline processes data from RNASeq geared towards transcription profiling, to generate a differential gene expression analysis.

### Major Steps:
1. Raw paired end read data is downloaded from the Genelab Data Repository
1. The raw reads undergo FastQC and MultiQC to allow the researcher to check raw read quality
1. The reads are trimmed by Trim-Galore. FastQC and MultiQC are repeated for the trimmed reads.
1. The Ensembl genome and gene annotations are downloaded
1. RSEM and STAR references are built using the Ensembl reference data
1. Trimmed reads are aligned to STAR reference genome
1. RSEM is used to count aligned reads per gene (and isoform)
1. An R script employs the library DESeq2 to tables of both normalized and raw counts.  Additionally, other related output is created including PCA and a statistical analysis of the counts data.

## Installation

This pipeline requires both Nextflow and Conda to be installed.

1. Ensure Nextflow is installed and available (tested with Nextflow version 20.10.0 build 5430):
<https://www.nextflow.io/docs/latest/getstarted.html>

1. Ensure Conda is installed and available (tested with Conda 4.8.3):
<https://docs.anaconda.com/anaconda/install/>

## Usage

#### Help Menu
```bash
nextflow run J-81/Nextflow_RCP -r test-awg-approved --help
```

#### Truncated Read and Genome Subsample Run
```bash
nextflow run J-81/Nextflow_RCP -r test-awg-approved --gldsAccession GLDS-194 --ensemblVersion=96 --truncateTo=300000 --genomeSubsample=19
```

### Optional Recommended Parameters

#### Nextflow Tower (Web-Based  Realitime Workflow Monitoring)
1. Register an account and generate a Nextflow tower token: [https://tower.nf](https://tower.nf)
1. Export token
```bash
export TOWER_ACCESS_TOKEN=<YOUR_TOKEN>
```
1. Add to your nextflow command
``` bash
nextflow run ... -with-tower
```

#### Resuming failed runs
1. Add to your nextflow command
``` bash
nextflow run ... -resume
```
**NOTE: For resume to correctly find cached successful tasks, you must run in the same directory as prior runs**


#### Reuse Shared References
- This saves resources by reusing previously downloaded reference files and builds.
1. Use shared references.
``` bash
nextflow run ... --storeDirPath='/data2/JO_Internship_2021/.References'
```


## License
[MIT](https://choosealicense.com/licenses/mit/)
