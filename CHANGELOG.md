# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).


## 0.0.4 - 2024-05-17

### Added
- support for single end data and long read data
- support for running no preprocessing of reads

## 0.0.3 - 2024-05-16

### Added
- `.tsv` based summaries of `breseq` output
- example `htcondor`-based running profile for `snakemake`>= 8.0

### Changed
- test module now creates data where `breseq` can detect deletions

## 0.0.2 - 2024-05-08

### Added
- De novo assembly with unicycler
- Comparison of de novo assembly to reference with quast

### Changed
- Updated pacakge versions to >= minimum requirements

## 0.0.1 - 2021-10-19

### Added
- A rule to rename the breseq `.gd` and `.vcf` output based on sample name

### Changed
- Updated `workflow/envs/multiqc_config.yaml` to better handle sample names
  across modules

### Bug fixes
- Fixed a bug in the `quality_control` rule where input samples could not be
  found
