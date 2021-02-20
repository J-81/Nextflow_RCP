# FORKED FROM J-81/masterProject

This repo includes the nextflow implemention of the RNASeq pipeline.  Note: the default output paths are designed to match the ones found in GLDS-194 as prepared on Oberyn.

---

############### Original README.md below ###############



# NASA Pipeline: GL-DPPD-7101-C

This is a nextflow implementation of a NASA bioinformatics pipeline for Jonathan Oribello's Bioinformatics Master's Project.  

This pipeline processes RNA_Seq, geared towards transcription profiling, to generate a differential gene expression analysis.

### Major Steps:
1. Raw paired end read data is downloaded from the Genelab Data Repository
1. The raw reads undergo FastQC and MultiQC to allow the researcher to check raw read quality
1. The reads are trimmed by Trim-Galore. FastQC and MultiQC are repeated for the trimmed reads.
1. The Mouse genome and gene annotations are downloaded from Ensembl: http://uswest.ensembl.org/Mus_musculus/Info/Index
1. RSEM and STAR references are built using the Ensembl data
1. Trimmed reads are aligned to STAR reference genome
1. RSEM is used to count aligned reads per gene (and isoform)
1. An R script employs the library DESeq2 to tables of both normalized and raw counts.  Additionally, other related output is created including PCA and a statistical analysis of the counts data.

## Installation

Ensure Conda is installed and available (tested with Conda 4.8.3):
<https://docs.anaconda.com/anaconda/install/>



Create conda environment using package main environment file and activate
```bash
conda env create -f envs/main.yml && conda activate main
```

## Usage

### Running Test Setup

** NOTE: AT THIS TIME, THIS WILL LIKELY FAIL ON LOCAL MACHINES AS BUILDING STAR REQUIRES >30 GB RAM **

Running with Tower monitoring (Optional, Highly Recommended):
- Login and setup here: [NextflowTower](https://tower.nf)
- **Ensure Nextflow environment variables are set**

```bash
nextflow pull J-81/masterProject && nextflow run J-81/masterProject -r dev -profile test  -with-tower
```

Running with without Tower monitoring (No Extra Setup):

```bash
nextflow pull J-81/masterProject && nextflow run J-81/masterProject -r dev -profile test
```

## Contributing
Outside contribution directly to this implementation is not welcome at this time.

## License
[MIT](https://choosealicense.com/licenses/mit/)
