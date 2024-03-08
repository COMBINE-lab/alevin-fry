# Changelog

Changelog for alevin-fry

## [0.9.0](https://github.com/COMBINE-lab/alevin-fry/compare/v0.8.2...v0.9.0) (2024-03-08)


### Features

* working with libradicl-0.8.2-pre ([12d088c](https://github.com/COMBINE-lab/alevin-fry/commit/12d088c836b4a75ac50b72aeb915888f900768d5))


### Bug Fixes

* MD to RST link formatting in overview ([672c987](https://github.com/COMBINE-lab/alevin-fry/commit/672c9879fc18fd9be9b2b6d08362bedeed6dd344))

## [0.8.2](https://github.com/COMBINE-lab/alevin-fry/compare/v0.8.1...v0.8.2) (2023-06-29)


### Bug Fixes

* update deps, make clippy happy ([ac1a316](https://github.com/COMBINE-lab/alevin-fry/commit/ac1a316001103f39a173b8a066c62810a8724875))

## [0.8.1](https://github.com/COMBINE-lab/alevin-fry/compare/v0.8.0...v0.8.1) (2023-01-12)


### Features

* **update_deps:** update dependencies ([a77c96e](https://github.com/COMBINE-lab/alevin-fry/commit/a77c96e162758e8cf5f4e509263216158bb580c9))


### Bug Fixes

* address [#95](https://github.com/COMBINE-lab/alevin-fry/issues/95) ([ad8ce7d](https://github.com/COMBINE-lab/alevin-fry/commit/ad8ce7d69a619d5f6d8bff8156483186293458e3))

## [0.8.0](https://github.com/COMBINE-lab/alevin-fry/compare/v0.7.0...v0.8.0) (2022-10-11)


### ⚠ BREAKING CHANGES

* fix parsing of force-cells and expect-cells

### Bug Fixes

* fix parsing of force-cells and expect-cells ([aa0f4ab](https://github.com/COMBINE-lab/alevin-fry/commit/aa0f4abb388e464491652fb2e1b682e33b1df05c))

## [0.7.0](https://github.com/COMBINE-lab/alevin-fry/compare/v0.6.0...v0.7.0) (2022-08-02)


### ⚠ BREAKING CHANGES

* try to force release-please

### Miscellaneous Chores

* try to force release-please ([4df099d](https://github.com/COMBINE-lab/alevin-fry/commit/4df099d9ee9e7c8cc81ef826d8723c7b6a5453ae))

## [0.6.0](https://github.com/COMBINE-lab/alevin-fry/compare/v0.5.1...v0.6.0) (2022-06-01)

### ⚠ BREAKING CHANGES

* **cmd_interface:** Add a (hidden) --umi-edit-dist option that take a value informing the underlying resolution algorithm of which edit distances to consider for collapse among potentially colliding UMIs. Right now, 0 works with all methods, while 1 only works with parsimony(-gene) and parimony(-gene)-em. The default remain as before, 1 for parsimony(-gene) and parsimony(-gene)-em, and 0 for all other methods. If the user attempts to set an unsupported edit distance for a method, the program will complain and exit.

### Features

* add hidden large PUG threshold cmd line option ([24fca6b](https://github.com/COMBINE-lab/alevin-fry/commit/24fca6b647a4686757a89f67e05808a151c6d231))
* add release-please ([ec4678b](https://github.com/COMBINE-lab/alevin-fry/commit/ec4678b7aa576daf1b41798d6d0614b09e08bbab))
* add toy run in github actions ([ea04c85](https://github.com/COMBINE-lab/alevin-fry/commit/ea04c855dbe556de2cdff324a1e02e6662656590))
* add toy run in github actions ([722d850](https://github.com/COMBINE-lab/alevin-fry/commit/722d850a70d6872510d9f5056e4f9a610f0b463b))
* cmdline file and directory validation ([9a21ee4](https://github.com/COMBINE-lab/alevin-fry/commit/9a21ee4c9ce0e06e63fee9f27f04fcf318a738b9))
* **resolution:** Add exact UMI only dedup mode to parsimony ([ac61b1d](https://github.com/COMBINE-lab/alevin-fry/commit/ac61b1d169db6b231f598c20348107f28e62d7dc))
* **resolution:** Add parsimony-gene and parsimony-gene-em modes ([1670faa](https://github.com/COMBINE-lab/alevin-fry/commit/1670faa4154042d2d855921a5a9df899c61ac5fa))
* **resolution:** add USA support to parsimony and parsimony-em ([8acc61e](https://github.com/COMBINE-lab/alevin-fry/commit/8acc61e53448ab326c2f61f3507c170b9448ccba))
* **resolution:** Improve parsimony in USA mode ([0f6eaa9](https://github.com/COMBINE-lab/alevin-fry/commit/0f6eaa940ab249cdf453a2868904e80a8b1d9383))
* toy run in github actions ([8d160cb](https://github.com/COMBINE-lab/alevin-fry/commit/8d160cbc76defc75da7a8073a853403ca848a7c1))
* toy run in github actions ([7331542](https://github.com/COMBINE-lab/alevin-fry/commit/7331542c4c71d53316d1781a783962307a596824))
* toy run in github actions ([f343bda](https://github.com/COMBINE-lab/alevin-fry/commit/f343bda9f121d1add7ff98cb7aa1deca7d0fd3b4))

### Bug Fixes

* **resolution:** Fix indexing in alternative resolution for usa parsimony ([8ebd3be](https://github.com/COMBINE-lab/alevin-fry/commit/8ebd3bed6caf74786d149899814c8715c994c041))
* **resolution:** revert prefer splicing heuristic ([4806146](https://github.com/COMBINE-lab/alevin-fry/commit/4806146394767bdbe2256ab8efa3c53a5f903c11))
* update compare_counts.py to user newer pyroe ([9a812a7](https://github.com/COMBINE-lab/alevin-fry/commit/9a812a7f8b57e42dce11a42983114311670856a4))

### Code Refactoring

* **cmd_interface:** remove --pug-exact-umi add --umi-edit-dist ([7145d64](https://github.com/COMBINE-lab/alevin-fry/commit/7145d64c2cabf8afd087dff2e4acb526d09a3bcb))

## [0.4.3] - 2021-11-11

### Fixed

- Fixed a bug that prevented the 1-edit rescue (in `generate-permit-list`) for barcodes of odd length, when using the unfiltered permit-list filtering mode. Thanks to @Gaura for helping to find and diagnose the issue.

## [0.4.2] - 2021-10-16

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

[unreleased]: https://github.com/COMBINE-lab/alevin-fry/compare/v0.4.3...HEAD
[0.4.3]: https://github.com/COMBINE-lab/alevin-fry/compare/v0.4.2...0.4.3
[0.4.2]: https://github.com/COMBINE-lab/alevin-fry/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/COMBINE-lab/alevin-fry/compare/v0.4.0...v0.4.1
