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
import os
from json.decoder import JSONDecodeError

# Overriding defaults with configuration from file
def get_config(override_name=None):
    """Get configuration parameters from ppn_config.py.

    If override_name is given, function updates non-default values."""
    params = defaults.copy()

    if override_name is not None:
        override_params = get_params(override_name)
        if not override_params:
            no_config_found_msg = ("Couldn't find overriding parameters in "
                                   "ppn_config. Key '{}' was unrecognised.").format(override_name)
            raise ValueError(no_config_found_msg)
        for key, group in override_params.items():
            try:
                params[key].update(group)
            except KeyError:
                params[key] = group
            except AttributeError:
                pass # Ignore keys that are not dictionaries

    ### Configuration input checks
    # For convenience
    data = params.get("data_source")
    runopt = params.get("run_options")
    ncopt = params.get("nowcast_options")

    _check_datasource(data)
    # Leadtimes generation
    _check_leadtime(params)

    # PGM-related checks
    if ("pgm" in data.get("importer", "").lower() and
            not params["output_options"].get("use_old_format", False)):
        raise ValueError(("Cannot write ODIM-compliant output when input data is in .pgm format! "
                         "Set 'output_options.use_old_format' to 'true'."))

    # Nowcast-method specific input checks
    if runopt.get("nowcast_method") in ["steps"]:
        # Timestep check
        d = data.get("timestep", None)
        n = ncopt.get("timestep", None)
        if n is None:
            ncopt["timestep"] = d
        elif n != d:
            error_msg_timestep_conflict = (
                "Conflicting options for 'timestep' in 'data_sources' and 'nowcast_options': "
                "{} and {}. If both are defined, they must be equal! (The parameter 'timestep' "
                "in nowcast_method is the time step for motion vectors, not for leadtimes.)"
            ).format(d, n)
            raise ValueError(error_msg_timestep_conflict)
        # kmperpixel check
        if ncopt.get("kmperpixel", None) is None and ncopt.get("vel_pert_method") in ["bps"]:
            raise ValueError("Configuration error: kmperpixel is required")

    # Expand ~ in paths, if any
    params["data_source"]["root_path"] = os.path.expanduser(params["data_source"]["root_path"])
    params["output_options"]["path"] = os.path.expanduser(params["output_options"]["path"])
    params["logging"]["log_folder"] = os.path.expanduser(params["logging"]["log_folder"])
    params["callback_options"]["tmp_folder"] = os.path.expanduser(params["callback_options"]["tmp_folder"])

    # Resolve relative paths, if any
    params["callback_options"]["tmp_folder"] = os.path.join(params["output_options"]["path"],
                                                            params["callback_options"]["tmp_folder"])

    return params

def _check_datasource(ds):
    """Check config file for errors in data_source definition"""
    if not isinstance(ds, dict):
        raise ValueError("Configuration error: mandatory option group 'data_source' is missing!")

    required_keys = ["root_path", "path_fmt", "fn_pattern", "fn_ext", "importer",
                     "timestep", "importer_kwargs"]  # importer_kwargs may be empty
    # Go through all keys and display all missing keys
    missing = []
    for key in required_keys:
        value = ds.get(key, None)
        if value is None:
            missing.append(key)
    if missing:
        raise ValueError("Following required parameters are missing from data_source group: {}".format(missing))
    # Type checks
    if not isinstance(ds.get("timestep"), int):
        raise TypeError("Configuration error in data_sources: timestep must be an integer")
    if not isinstance(ds.get("importer_kwargs"), dict):
        raise TypeError('Configuration error in data_sources: importer_kwargs must be an object '
                        '({"key": value}). It can be empty.')

# FIXME: Logic could be simplified
def _check_leadtime(params):
    # For convenience
    runopt = params.get("run_options")
    leadtimes = runopt.get("leadtimes")
    max_leadtime = runopt.get("max_leadtime")
    nowcast_timestep = runopt.get("nowcast_timestep")
    input_timestep = params.get("data_source", dict()).get("timestep")

    # leadtimes is not given, but timestep and maximum is
    if leadtimes is None:
        if None in [max_leadtime, nowcast_timestep]:
            raise ValueError("Need to set both 'max_leadtime' and 'nowcast_timestep' if 'leadtimes'"
                             " is not given")
        leadtimes = int(max_leadtime / nowcast_timestep)
        runopt["leadtimes"] = leadtimes

    # leadtimes is given as number ("this many leadtimes"), check the input data timestep
    if isinstance(leadtimes, int):
        # if input data timestep == nowcast timestep or nowcast timestep is None, pass
        # else generate a list
        if nowcast_timestep is not None and nowcast_timestep != input_timestep:
            runopt["leadtimes"] = [nowcast_timestep + i*nowcast_timestep for i in range(leadtimes)]

    # if leadtimes is a list, pass (might need to include a check for valid values within list)
    elif isinstance(leadtimes, (list, tuple)):
        pass
    else:
        raise ValueError("Cannot figure out leadtimes")
    # else raise error?

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
    name = os.path.splitext(name)[0] + ".json"
    if os.path.exists(name):
        cfname = name
    else:
        cfpath = os.path.realpath(__file__)
        cfpath = os.path.split(cfpath)[0]
        cfname = os.path.join(cfpath, "config/", name)
    print(f"Using PPN configuration from {cfname}")
    try:
        with open(cfname, "r") as f:
            params = json.load(f)
    except FileNotFoundError as exc:
        file_missing_msg = ("Cannot find requested config '{}'! Are the filename and path correct?"
                            "\n{}").format(os.path.splitext(name)[0], cfname)
        raise OSError(file_missing_msg) from exc
    except JSONDecodeError as exc:
        raise RuntimeError("Could not decode config file. Is it valid JSON file?") from exc

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

def dump_empty():
    """Utility function for generating an empty configuration file."""
    empty = {key: {} for key in defaults.keys()}
    dump_params_to_json(empty, "empty")

# Default parameters for PPN, other dictionaries should override these parameters
# using dict.update() method.
defaults = {
    # Option groups in alphabetical order
    "data_options": {
        "zr_a": 223,
        "zr_b": 1.53,
        "rain_threshold": 8,  # In data units
    },

    "data_source": {
    },

    "logging": {
        "write_log": False,
        "log_level": logging.INFO,
        "log_folder": "/tmp",
    },

    "motion_options": {
    },

    "nowcast_options": {
        # Default to the nowcast method defaults
        # n_ens_members = 24,
        "n_cascade_levels": 6,  # pysteps default for steps
        "fft_method": "pyfftw",
        "vel_pert_method": "bps",  # pysteps default, requires kmperpixel to be set
        "vel_pert_kwargs": {
            # lucaskanade/fmi values given in pysteps.nowcasts.steps.forecast() method documentation
            "p_par": [2.20837526, 0.33887032, -2.48995355],
            "p_perp": [2.21722634, 0.32359621, -2.57402761],
        },
        "domain": "spectral",  # pysteps default

        "num_workers": 6, # pysteps defaults to 1, we want more.
        "seed": None,  # pysteps default

        # Previously hardcoded values
        "extrap_method": "semilagrangian",
        "noise_method": "nonparametric",
        "ar_order": 2,
        "mask_method": "incremental",
        # Required parameters that need to be calculated or defined
        #"kmperpixel": # required for motion perturbation method (bps)
        #"timestep": # "Timestep" of MOTION VECTOR, get from input data
        #"R_thr":  # must be in decibel units...
    },

    "output_options": {
        #
        "path": "/tmp",
        #
        "store_motion": True,
        "store_perturbed_motion": False,
        "store_ensemble": True,
        "store_deterministic": True,
        #
        "as_quantity": None, # None == same as input. Other valids: DBZH, RATE (see ODIM standard)
        "convert_to_dtype": None,  # None == same as input. Otherwise provide a valid string for numpy.dtype()
        #"scaler": 10, # Deprecated
        #"scale_zero": "auto", # Deprecated
        "gain": 0.1,
        "offset": -32,  # "auto" for data minimum
        # "input" values are encoded, numbers given here are not
        "set_undetect_value_to": "input", # a number or "input". input == read from input (units converted if needed)
        "set_nodata_value_to": "default",  # A number, "default", "max_int", "min_int", or "nan" (floats only).
        #
        "write_leadtimes_separately": False, # Store each leadtime after calculating it instead of everything at the end
        "write_asap": True,
        "use_old_format": False,  # Remove when postprocessing can use ODIM format
    },

    "run_options": {
        "leadtimes": 12,  # int = number of timesteps, list of floats = forecast for these lead times
        # if leadtimes is not a list and nowcast_timestep != input timestep, use these to make it into one
        "nowcast_timestep": None, # optional, default to input timestep
        "max_leadtime": 60, # optional, used only if "leadtimes" is None
        # What is calculated
        "run_deterministic": True,
        "run_ensemble": True,
        "regenerate_perturbed_motion": False,  # Re-calculate the perturbed motion fields used for pysteps nowcasting
        #
        "num_prev_observations": 3,
        # Methods
        "motion_method": "lucaskanade",
        "nowcast_method": "steps",
        "deterministic_method": "extrapolation",

        #
        "forecast_as_quantity": "DBZH",  # Input data is converted to this before nowcasting
        "steps_set_no_rain_to_value": -10,  # In forecast quantity units
    },

    # Used when writing ensemble nowcasts after each timestep with callback function
    "callback_options": {
        "tmp_folder": "tmp",  # relative to output_options.path (or absolute path)
    }
}

# Test cases
if __name__ == '__main__':
    from pprint import pprint
    # Test custom parameter updating
    print("Updating default parameters with custom config")
    cfg = get_config("test/customised")
    pprint(cfg)
