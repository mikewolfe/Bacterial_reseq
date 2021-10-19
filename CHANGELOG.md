# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## 0.0.1 - 2021-10-19

### Added
- A rule to rename the breseq `.gd` and `.vcf` output based on sample name

### Changed
- Updated `workflow/envs/multiqc_config.yaml` to better handle sample names
  across modules

### Bug fixes
- Fixed a bug in the `quality_control` rule where input samples could not be
  found
