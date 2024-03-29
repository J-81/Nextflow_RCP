# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [rc1.0.5] - 2022-06-02

### Added

#### STAR Alignment

- V&V implemented akin to RSEM counts table
  - Ensure counts correctly aggregated from {sample}_ReadsPerGene.out.tab

### Changed

#### RSEM Index naming

- Uses {Oorg}.grp for file prefix (e.g. Mmus.grp / Mmus_ERCC.grp)

#### RSEM Count

- Unnormalized_Counts.csv renamed to RSEM_Unnormalized_Counts.csv to disambiguate compared to STAR_Unnormalized_Counts.csv

#### Staging

- raw fastq.gz files will now be named as they were originally named in the GLDS repo
  - Files derived from those original raw fastq.gz will be renamed to a standard format (as opposed to renaming the original raw files in a standard format)

### Fixed

#### DESEQ2 Annotations

- Modified to use database files generated in-house.
  - This addresses a bug in STRINGdb and may be reverted once the STRINGdb bug fix is released

#### Post Processing

- md5sum table is not synced with assay table
  - Now correctly synced

## [rc1.0.4] - 2022-06-02

### STAR Alignment

#### Added

- New published files
  - Output directory: 02-STAR_Alignment
    - STAR_NumNonZeroGenes.csv
      - A dataset wide tabulation of the number of genes with non zero counts.
    - STAR_Unnormalized_Counts.csv
      - A dataset wide tabulation of gene counts output from STAR.

#### Fixed

- STAR_Unnormalized_Counts.csv file is published as originally intended

## [rc1.0.3] - 2022-04-30

### Added

#### STAR Alignment

- Will now output quantification 'GeneCounts' (in addition to TranscriptomeSAM)

### Changed

### Fixed

### Removed

## [0.1.1-beta] - 2022-03-09
### Fixed
  - Python scripts set to executable

## [0.1.0-beta] - 2022-01-21

[0.1.1-beta]: https://github.com/J-81/Nextflow_RCP/compare/0.1.1-beta...0.1.0-beta
