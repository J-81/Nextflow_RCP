# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Planned
### Added
#### STAR Alignment
  - STAR_counts_table.R process that summarizes new GeneCounts output in output file: STAR_Unnormalized_Counts.csv

### Changed
#### RSEM Count
  - Unnormalized_Counts.csv renamed to RSEM_Unnormalized_Counts.csv to disambiguate compared to STAR_Unnormalized_Counts.csv

### Changed
#### Staging
  - raw fastq.gz files will now be named as they were originally named in the GLDS repo
    - Files derived from those original raw fastq.gz will be renamed to a standard format (as opposed to renaming the original raw files in a standard format)

## [Unreleased] Target Release: 0.1.2-beta
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
