# Changelog

Changelog for alevin-fry

## [Unreleased]

Currently no unreleased features

## [0.4.2] - in progress

### Added

- Support for USA mode to the `infer` command by passing the `--usa` flag.
### Changed

- Large _internal_ code re-organization, moving most alevin-fry related functionality out of the `libradicl` crate / library.
- Some details of how `cr-like-em` works in USA mode. Instead of each splicing status for each gene being completely independent, the expecation of assignment for spliced and unspliced also depend on ambiguous abundance, and ambiguous abundance depends on both spliced and unspliced.
### Fixed
- [Issue #25](https://github.com/COMBINE-lab/alevin-fry/issues/25) where cells with only (and too few) highly-multimapping reads could sometimes prevent quantification from completing successfully when using the `cr-like-em` mode.  Now, such cells are instead flagged 
and their identifiers are output in `quant.json`.  Generally, these cells will have no expressed genes in the corresponding count matrix.
## [0.4.1] - 2021-07-22

This is a minor release, intended mostly to bump some version dependencies and to address [Issue #22](https://github.com/COMBINE-lab/alevin-fry/issues/22).

### Changed

- Changed the name of the JSON file written by the `quant` command from `meta_info.json` to `quant.json` to match other commands
- Updated versions of crates in dependencies for libradicl and alevin-fry

### Removed

- The `quant` command no longer writes a `cmd_info.json` file

[unreleased]: https://github.com/COMBINE-lab/alevin-fry/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/COMBINE-lab/alevin-fry/compare/v0.4.0...v0.4.1
