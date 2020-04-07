# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project DOES NOT adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Instead, version numbering will be a running number. Each major number increases when a dependency module is upgraded, added, or removed (e.g. pysteps). Each minor number increases when a research version is updated, either by RAS or PAK.

## [v3.0] - 2020-01-29
### Added
- Folder for configuration files
- Utility functions for generating a configuration json file
- Parameters `NORAIN_VALUE` and `RAIN_THRESHOLD` for thresholding
- MIT License

### Changed
- Configurations are now read from json files
- Updated default configuration based on Pysteps article (Pulkkinen et al. 2019)
- Updated pysteps version to 1.1.1
- `LOG_FOLDER`, `LOG_LEVEL` and `OUTPUT_PATH` configuration parameters are no longer included in output metadata

### Removed
- Configuration parameters `DBZ_MIN`, `DBZ_THRESHOLD`, `R_MIN`, `R_THRESHOLD` (Replaced with new parameters).

### Fixed
- Program tried to write a log file before logging was initialised, when program couldn't find a given configuration. Now program raises an error instead.
- Input data is now thresholded properly
- Now it is possible to output values in mm/h, even if calculations are done in dBZ

## [v2.1] - 2019-12-20
### Added
- New parametrisations and configurations for operational use
- Usage and Configuration descriptions to readme

### Changed
- Now `pystepsrc` file is tracked, previously it was ignored
- Renamed "dev" parametrisation to "test"
- Updated default parameters: now generates 25 ensemble members for up to 120 minutes
- Updated docstrings
- Updated .gitignore
- Improved markdown syntax usage in changelog and readme

### Removed
- Parametrisations used in research

### Fixed
- Dangerous default argument value in `utils.store_timeseries()`
- Removed unused import from `utils`
- Version tag links in changelog

## [v2.0] - 2019-09-02
### Added
- Changelog
- New dependency: pysteps (1.0.0)
- Now it is possible to set motion perturbation parameters
- Default motion perturbation parameters for Lucas-Kanade method
- STEPS nowcasts can now be calculated in dBZ, no need to convert to rain rate first
- Functionality for regenerating perturbed motion fields (for a given `SEED`)
- Projection metadata is now stored in output file
- Configuration parameters and their values are now stored in output file
- Added new configuration parameters for controlling program flow
- Added configuration for heavy rainfall event in Helsinki on 23 Aug 2019

### Changed
- Renamed readme.md to README.md
- Deterministic nowcast uses now extrapolation
- Old deterministic nowcast is renamed as "unperturbed nowcast"
- `utcnow_floored()` function was moved to `utils` module
- Edited "dev" and "verification" parameters a bit

### Removed
- Extra dependencies that were unused or no longer needed

## [v1.1] - 2019-05-06
### Added
- New parameters: `SEED`, `ZR_A`, `ZR_B`, `SCALE_ZERO`, `SCALER`, `FIELD_VALUES`
- Now it is possible to set random number generator seed via `SEED` parameter
- Now R(Z)-relation coefficients _a_ and _b_ can be set via parameters `ZR_A` and `ZR_B`
- Nowcasts can be stored either in dBZ or mm/h units. This is chosen via `FIELD_VALUES` parameter
- Scaling parameter in output file can now be set via parameter `SCALER`
- Scaling offset in output file can now be set via `SCALE_ZERO` parameter
- Program uses nowcast data array's minimum non-NaN value as scaling offset, unless given via configuration parameters
- Parameters `SEED` and `FIELD_VALUES` are stored in output file`s metadata

### Changed
- "Valid for" attribute in output file is now an integer (was: string)
- `OUTPUT_TIME_FORMAT` parameter no longer affects the "Valid for" attribute

## [v1.0] - 2019-01-22
- First prototype version given to operational testing

[v3.0]: https://github.com/fmidev/fmippn/compare/v2.1...v3.0
[v2.1]: https://github.com/fmidev/fmippn/compare/v2.0...v2.1
[v2.0]: https://github.com/fmidev/fmippn/compare/v1.1...v2.0
[v1.1]: https://github.com/fmidev/fmippn/compare/v1.0...v1.1
[v1.0]: https://github.com/fmidev/fmippn/releases/tag/v1.0

