"""Convert configuration files to the new format"""
import os
import sys
import json
import shutil
from pprint import pprint

import pysteps
from pysteps import rcparams as pystepsrc

def mapping(old_config):
    """Map old config parameters to new ones with default parameters from old style."""

    vel_pert_defaults = {
        # lucaskanade/fmi values given in pysteps.nowcasts.steps.forecast() method documentation
        "p_par": [2.20837526, 0.33887032, -2.48995355],
        "p_perp": [2.21722634, 0.32359621, -2.57402761],
    }

    old_to_odim = {
        'dbz': 'DBZH',
        'rrate': 'RATE',
    }

    dataopt = {
        "zr_a": old_config.get("ZR_A", 223.),
        "zr_b": old_config.get("ZR_B", 1.53),
        "rain_threshold": old_config.get("RAIN_THRESHOLD"), # in data units, e.g. 8 dBZ or 0.1 mm/h
    }
    outopt = {
        #"path": old_config.get("OUTPUT_PATH", '/tmp'), # Default to new style
        "as_quantity": old_to_odim[old_config.get("FIELD_VALUES", "dbz")], # Store values as rrate or dbz

        "scaler": old_config.get("SCALER", 100),
        "scale_zero": old_config.get("SCALE_ZERO", "auto"),
        "store_ensemble": old_config.get("STORE_ENSEMBLE", True), # Write each ensemble member to output
        "store_deterministic": old_config.get("STORE_DETERMINISTIC", True), # Write det_fct to output
        "store_motion": old_config.get("STORE_MOTION", True), # Write deterministic motion to output
        "store_perturbed_motion": old_config.get("STORE_PERTURBED_MOTION", True), # Write motion for each ensemble member to output

        "use_old_format": True,
    }
    # This cannot be None, but it was old default
    old_path = old_config.get("OUTPUT_PATH")
    outopt["path"] = old_path if old_path is not None else '/tmp'

    logopt = {
        "write_log": old_config.get("WRITE_LOG", False),
        "log_level": old_config.get("LOG_LEVEL", 20),  # see logging module for values
        "log_folder": old_config.get("LOG_FOLDER", "../logs"),
    }

    runopt = {
        "leadtimes": old_config.get("NUM_TIMESTEPS", None), # or MAX_LEADTIME/NOWCAST_TIMESTEP
        "nowcast_timestep": old_config.get("NOWCAST_TIMESTEP", 5), # optional, default to input timestep
        "num_prev_observations": old_config.get("NUM_PREV_OBSERVATIONS", 3),
        "max_leadtime": old_config.get("MAX_LEADTIME", 120),

        "run_deterministic": old_config.get("GENERATE_DETERMINISTIC", True),
        "run_ensemble": old_config.get("GENERATE_ENSEMBLE", True),
        "regenerate_perturbed_motion": old_config.get("REGENERATE_PERTURBED_MOTION", False),  # Re-calculate the perturbed motion fields used for pysteps nowcasting

        "motion_method": old_config.get("OPTFLOW_METHOD", "lucaskanade"),
        "nowcast_method": "steps",  # Hardcoded in PPN
        "deterministic_method": "extrapolation", # Hardcoded in PPN

        "forecast_as_quantity": old_to_odim[old_config.get("VALUE_DOMAIN", "dbz")],
        "steps_set_no_rain_to_value": old_config.get("NORAIN_VALUE"),  # In forecast quantity units

    }

    # Unwanted options should be removed to default to nowcast method's values
    ncopt =  {
        "kmperpixel": old_config.get("KMPERPIXEL"),

        "timestep": old_config.get("NOWCAST_TIMESTEP", 5),  # "Timestep" of MOTION VECTOR
        "fft_method": old_config.get("FFT_METHOD", "pyfftw"),
        "n_ens_members": old_config.get("ENSEMBLE_SIZE", 24),
        "n_cascade_levels": old_config.get("NUM_CASCADES", 8),
        "num_workers": old_config.get("NUM_WORKERS", 6),
        "domain": old_config.get("CALCULATION_DOMAIN", "spectral"),  # default used to be "spatial"
        "vel_pert_method": old_config.get("VEL_PERT_METHOD", "bps"),
        "vel_pert_kwargs": old_config.get("VEL_PERT_KWARGS", vel_pert_defaults),
        "seed": old_config.get("SEED", None),

        # Previously hardcoded values
        #"extrap_method": "semilagrangian",
        #"noise_method": "nonparametric",
        #"ar_order": 2,
        #"mask_method": "incremental",
    }

    new_config = {
        # Used to be in pystepsrc file, now included in this config
        "data_source": pystepsrc["data_sources"].get(old_config["DOMAIN"], dict()),
        # Passed to optical flow method. Currently no options are set.
        "motion_options": {},

        # Input data related
        "data_options": dataopt,
        # Output data related
        "output_options": outopt,
        # Used methods etc.
        "run_options": runopt,
        # Passed to nowcast method
        "nowcast_options": ncopt,
        # Finally, logging stuff
        "logging": logopt,
    }

    return new_config

if __name__ == '__main__':
    cwd = os.getcwd()
    rcpath = os.path.dirname(pysteps.config_fname())
    if cwd != rcpath:
      raise RuntimeError("Please run this script in the same folder as your pystepsrc!\n"
                         "    E.g. 'python tools/convert_config.py test.json'")
    if len(sys.argv) < 2:
        raise ValueError("Please give config file name (e.g. test.json)")
    config_folder = "./config"
    backup_folder = "./tools/bak"
    print(f"Creating a backup in folder {backup_folder}")
    os.makedirs(backup_folder, exist_ok=True)
    cfgname = os.path.join(config_folder, sys.argv[1])
    shutil.copy2(cfgname, backup_folder)

    print(f"Trying to convert {cfgname}")
    with open(cfgname, 'r') as f:
        old = json.load(f)
    print(f'Reading data_sources from pystepsrc for DOMAIN: "{old["DOMAIN"]}"')
    if pystepsrc["data_sources"].get(old["DOMAIN"]) is None:
        print("Could not find proper data_sources from pystepsrc, please fill it manually to "
              "converted file.")
    new = mapping(old)
    # ~ pprint(old)
    # ~ pprint(new)
    # ~ pprint(json.dumps(new, indent=2))
    print(f"Overwriting old file {cfgname}")
    with open(cfgname, 'w') as out:
      json.dump(new, out, indent=2)
    print("Done.")
