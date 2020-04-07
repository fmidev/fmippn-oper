"""Parameter configuration for FMI-PPN

Dictionary `defaults` contains default parameter settings for FMI-PPN.

Other dictionaries should contain non-default parameters. Program should
update defaults with non-default dictionary.

Utility method `get_params(name)` returns a dictionary.

Adding new parametrizations:
  1. Create a new json file in config folder
  2. Add new parameter values in a json dictionary object

TIP: To easily create a json file, dump the defaults from this file and modify
the resulting json file. Just remember to rename the json file! Easy way of
generating the file is to run following command on command line:

/path/to/fmippn/source$ python -c "import ppn_config; ppn_config.dump_defaults()"
"""
import logging
import json

def get_params(name):
    """Utility function for easier access to wanted parametrizations.

    Two parametrisations should be provided: `defaults` and `test`.
    `defaults` contains all possible parameters and default values for them.
    The default values are used unless explicitly overriden.

    `test` contains parametrisations for short nowcast that generates all
    different nowcasts. Suitable for quickly testing different parameters and
    for development purposes.

    Input:
        name -- name of the parametrization dictionary

    Return dictionary of parameters for overriding defaults. Non-existing
    names will return an empty dictionary.
    """
    if name == "defaults":
        params = defaults
    else:
        try:
            with open(f"config/{name}.json", "r") as f:
                params = json.load(f)
        except FileNotFoundError:
            params = dict()

    return params

def dump_params_to_json(params, config_name):
    """Utility function for dumping the used parametrisations to a new config file.

    Input:
        params -- a Python dictionary containing the parametrisations
        config_name -- name for the new configuration file

    The parametrisations will be written to ./config/{config_name}.json
    """
    with open(f"config/{config_name}.json", "w") as f:
        json.dump(params, f, indent=2)

def dump_defaults():
    """Utility function for generating a configuration file with default values.
    Useful for creating new configurations using Python shell."""
    dump_params_to_json(defaults, "defaults")

# Default parameters for PPN, other dictionaries should override these parameters
# using dict.update() method.
defaults = {
    # Method selections
    "DOMAIN": "fmi", # See pystepsrc for valid data sources
    "OPTFLOW_METHOD": "lucaskanade",
    "FFT_METHOD": "pyfftw",
    "GENERATE_DETERMINISTIC": True,
    "GENERATE_ENSEMBLE": True,
    "GENERATE_UNPERTURBED": False,
    "REGENERATE_PERTURBED_MOTION": False,  # Re-calculate the perturbed motion fields used for pysteps nowcasting
    "VALUE_DOMAIN": "dbz",  # dbz or rrate
    # Z-R conversion parameters
    "ZR_A": 223.,
    "ZR_B": 1.53,
    # Nowcasting parameters
    "NUM_PREV_OBSERVATIONS": 3,
    "NOWCAST_TIMESTEP": 5,
    "MAX_LEADTIME": 120,
    "NUM_TIMESTEPS": None,
    "ENSEMBLE_SIZE": 24,
    "NUM_CASCADES": 8,
    "NUM_WORKERS": 6,
    "RAIN_THRESHOLD": 6.5,  # Roughly 0.1 mm/h using Z = 223 * R ** 1.53
    "NORAIN_VALUE": 1.5,  # Threshold minus 5 (dB) units. Roughly 0.04 mm/h using above
    "KMPERPIXEL": 1.0,
    "SEED": None,  # Default value in pysteps is None
    # Motion perturbation parameters
    # Set to VEL_PERT_KWARGS to `None` to use pysteps's default values
    "VEL_PERT_METHOD": "bps",
    "VEL_PERT_KWARGS": {
        # lucaskanade/fmi values given in pysteps.nowcasts.steps.forecast() method documentation
        "p_par": [2.20837526, 0.33887032, -2.48995355],
        "p_perp": [2.21722634, 0.32359621, -2.57402761],
    },
    # Storing parameters
    "FIELD_VALUES": "dbz",  # Store values as rrate or dbz
    "OUTPUT_PATH": None,  # None uses pystepsrc output path
    "OUTPUT_TIME_FORMAT": "%Y-%m-%d %H:%M:%S",
    "STORE_ENSEMBLE": True, # Write each ensemble member to output
    "STORE_UNPERTURBED": True,
    "STORE_DETERMINISTIC": True,  # Write det_fct to output
    "STORE_MOTION": True, # Write deterministic motion to output
    "STORE_PERTURBED_MOTION": True,  # Write motion for each ensemble member to output
    "SCALER": 100,
    "SCALE_ZERO": "auto",  # Value for "0" in scaled units. Set to "auto" or None for minimum value found before scaling
    # Logging parameters
    "WRITE_LOG": False,
    "LOG_LEVEL": logging.INFO,  # see logging module documentation for valid levels
    "LOG_FOLDER": "../logs",
}
