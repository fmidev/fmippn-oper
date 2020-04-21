# Finnish Meteorological Institute Probabilistic Precipitation Nowcasting system (FMI-PPN)
FMI-PPN is a modular weather-radar-based nowcasting system built for research and operational usage.

## Subsystems
### Precipitation motion
Current implementation of FMI-PPN uses Lucas-Kanade optical flow method to estimate precipitation movement from radar measurements.

### Nowcasting
Currently FMI-PPN uses [pysteps](https://pysteps.github.io) to generate ensemble nowcasts.

## Known issues and limitations
- Parameter `SEED` must be `None` or an integer between `0` and `2**32 - 1`. This is a limitation in `numpy.random`.
- OpenMPI conflicts with dask when both are installed, leading to significant decrease FMI-PPN performance.
  - Workaround is to set OpenMPI use only one thread. In Linux you can set environment variable `OMP_NUM_THREADS=1`.
  - Alternatively: uninstall `dask` from conda environment
- `pyfftw <0.12.0` is incompatible with `scipy 1.4`.

## Usage
### Installation
1. Install [conda](https://conda.io/en/latest/) (Miniconda is recommended)
2. Clone this repository
3. Change directory to FMI-PPN folder
4. Setup conda environment: `$ conda env create -f environment.yml`
5. Modify pystepsrc file: Add your data source configuration under `data_sources` (see pysteps documentation for details)
6. Replace the value for `DOMAIN` parameter under `defaults` with your data source configuration's name
7. Activate conda environment: `$ conda activate fmippn`
8. Run FMI-PPN with default settings: `$ python run_ppn.py`

### Running FMI-PPN
Before running FMI-PPN, you should configure pysteps (via `pystepsrc` file) and PPN by adding your parametrisations to a config file (see below).

How to run:
1. Activate your conda environment
2. Run FMI-PPN with your settings: `$ python run_ppn.py -c your_config_parametrisation`

## Configuration
### Adding new parametrisations
Create a `json` file with parameters you want to change from defaults and put it in `config` folder. The file's name (without file type extension) will be used to select the settings.

The `ppn_config.py` module has a utility function `dump_defaults()` for creating a configuration file based on default settings.

### Parametrisations
Parameter|Explanation|Default value
----|----|----
`CALCULATION_DOMAIN`|Choose if the nowcast calculation is performed in spatial (`spatial`) or spectral (`spectral`) domain. See also [pysteps documentation](https://pysteps.readthedocs.io/en/latest/generated/pysteps.nowcasts.steps.forecast.html#pysteps.nowcasts.steps.forecast).|`spectral`
`DOMAIN`|Data source used from pystepsrc|`fmi`
`ENSEMBLE_SIZE`|Number of ensemble members|`24`
`FFT_METHOD`|FFT method used in pysteps calculations|`pyfftw`
`FIELD_VALUES`|Select the units to store the nowcast (before scaling). Valid units are `dbz` for dBZ and `rrate` for mm/h.|`dbz`
`GENERATE_DETERMINISTIC`|Calculate extrapolation-only nowcast|`True`
`GENERATE_ENSEMBLE`|If `True`, then calculate ensemble members|`True`
`GENERATE_UNPERTURBED`|Calculate nowcast using pysteps, but without noise|`False`
`KMPERPIXEL`|Pixel size in kilometers|`1.0`
`LOG_FOLDER`|Path where log files should be stored|`../logs`
`LOG_LEVEL`|Logging level used by [Python Logging Module](https://docs.python.org/3/library/logging.html#levels)|`20`
`MAX_LEADTIME`|How long your nowcast will be (in minutes)|`120`
`NORAIN_VALUE`|Value assigned to dry pixels during thresholding. Units depend on `VALUE_DOMAIN` parameter (dBZ or mm/h). Must be less than `RAIN_THRESHOLD` value!|`1.5`
`NOWCAST_TIMESTEP`|Timestep between consecutive nowcast images|`5`
`NUM_CASCADES`|How many cascade levels are used in cascade decomposition by pysteps|`8`
`NUM_PREV_OBSERVATIONS`|Number of previous observations used in optical flow calculation|`3`
`NUM_TIMESTEPS`|How many timesteps will be calculated in nowcasts (If this setting is `None`, this value is automatically calculated based on `MAX_LEADTIME` and `NOWCAST_TIMESTEP` parameters) |`None`
`NUM_WORKERS`|Number of worker threads used in parallel computing|`6`
`OPTFLOW_METHOD`|Optical flow method (see pysteps docs)|`lucaskanade`
`OUTPUT_PATH`|If not `None`, then use this path to store output instead of setting in pysteprc|`None`
`OUTPUT_TIME_FORMAT`|Python datetime format for showing timestamps|`%Y-%m-%d %H:%M:%S`
`RAIN_THRESHOLD`|Thresholding value for rain. Pixels with values under this parameter are regarded as dry pixels. Units depend on `VALUE_DOMAIN` parameter (dBZ or mm/h).|`6.5`
`REGENERATE_PERTURBED_MOTION`|Re-calculate motion for each ensemble member (requires `SEED != None`) |`False`
`SCALER`|Scaling coefficient for output|`100`
`SCALE_ZERO`|Value for "0" after scaling the output. Setting to `"auto"` or `None` uses the minimum value found in data before scaling.|`auto`
`SEED`|Seed parameter for random number generation. Use `None` for unseeded nowcasts.|`None`
`STORE_DETERMINISTIC`|Write "deterministic" nowcast in output file|`True`
`STORE_ENSEMBLE`|Write all ensemble members in output file|`True`
`STORE_MOTION`|Write optical flow motion field in output file|`True`
`STORE_PERTURBED_MOTION`|Write optical flow motion for each ensemble member in output file|`True`
`STORE_UNPERTURBED`|Write "unperturbed" nowcast in output file|`True`
`VALUE_DOMAIN`|Choose if nowcasting is performed for data in dBZ or mm/h units. Valid parameters are `dbz` for dBZ and `rrate` for mm/h.|`dbz`
`VEL_PERT_KWARGS`|Parameters for velocity perturbation (see pysteps docs). Set this parameter to `None` to use default pysteps values.|`{'p_par': [2.20837526, 0.33887032, -2.48995355], 'p_perp': [2.21722634, 0.32359621, -2.57402761]}`
`VEL_PERT_METHOD`|Velocity perturbation method used in pysteps|`bps`
`WRITE_LOG`|If `True`, then generate a log file|`False`
`ZR_A`|Value for coefficient _a_ in R(Z) relation|`223.0`
`ZR_B`|Value for coefficient _b_ in R(Z) relation|`1.53`
