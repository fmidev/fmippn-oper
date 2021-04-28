"""
This is a version of FMI-PPN used for probability of lightning
nowcasts.
"""
import datetime as dt
import os
import random

import numpy as np
import h5py
import pysteps
from pysteps import rcparams as pystepsrc

import ppn_logger
import ppn_config
import utils

# Global object for storing and accessing configuration parameters
PD = dict()

def run(timestamp=None, config=None, **kwargs):
    """Main function for FMI-PPN.
    Input:
        timestamp -- timestamp of form YYYYMMDDHHMM (str)
                     If None, use latest composite (default=None)
        config -- Configuration parameter. If None, use defaults. (default=None)

    Optional keyword arguments:
        (none)
    """
    nc_fname = None

    PD.update(ppn_config.get_config(config))

    # Default is generate nowcasts from previous radar composite
    if timestamp is not None:
        startdate = dt.datetime.strptime(timestamp, "%Y%m%d%H%M")
    else:
        startdate = utils.utcnow_floored(increment=5)

    # Paths, importers etc.
    datasource = PD.get("data_source")
    # NOTE: This is for backwards compability, can be removed at some point
    if datasource is None:
        datasource = pystepsrc["data_sources"][PD["DOMAIN"]]

    # Used methods
    optflow = pysteps.motion.get_method(PD["run_options"]["motion_method"], **kwargs)
    nowcaster = pysteps.nowcasts.get_method(PD["run_options"]["nowcast_method"], **kwargs)
    deterministic_nowcaster = pysteps.nowcasts.get_method(PD["run_options"]["deterministic_method"], **kwargs)

    time_at_start = dt.datetime.today()

    run_options = PD["run_options"]

    input_files = get_filelist(startdate, datasource)

    #Korvaa tama Partion pygrib-koodilla!
    observations, obs_metadata = read_observations(input_files, datasource, importer)

    motion_field = optflow(observations, **PD.get("motion_options", dict()))

    deterministic, det_meta = generate_deterministic(observations[-1],
                                                         motion_field,
                                                         deterministic_nowcaster,
                                                         metadata=obs_metadata)

    #Output kirjoitus tahan
    

def get_filelist(startdate, datasource):
    """Get a list of input file names"""
    try:
        filelist = pysteps.io.find_by_date(startdate,
                                           datasource["root_path"],
                                           datasource["path_fmt"],
                                           datasource["fn_pattern"],
                                           datasource["fn_ext"],
                                           datasource["timestep"],
                                           num_prev_files=PD["run_options"]["num_prev_observations"])
    except OSError as pysteps_error:
        error_msg = "Failed to read input data!"
        log("error", f"OSError was raised: {error_msg}")
        # Re-raise so traceback is shown in stdout and program stops
        raise OSError(error_msg) from pysteps_error
    return filelist


def generate_deterministic(observations, motion_field, nowcaster, nowcast_kwargs=None,
                           metadata=None):
    """Generate a deterministic nowcast using semilagrangian extrapolation"""
    # Extrapolation scheme doesn't use the same nowcast_kwargs as steps
    if nowcast_kwargs is None:
        nowcast_kwargs = dict()
    forecast, meta = generate(observations, motion_field, nowcaster, nowcast_kwargs,
                              metadata)
    return forecast, meta



def prepare_data_for_writing(forecast):
    """Convert and scale ensemble and deterministic forecast data to uint16 type"""
    if forecast is None:
        return None, dict()

    # Store data in integer format to save space (float64 -> uint16)
    store_dtype = 'uint16'
    store_nodata_value = np.iinfo(store_dtype).max if store_dtype.startswith('u') else -1
    scaler = PD["output_options"]["scaler"]
    scale_zero = PD["output_options"]["scale_zero"]
    if scale_zero in [None, "auto"]:
        scale_zero = np.nanmin(forecast)
    prepared_forecast = utils.prepare_fct_for_saving(forecast, scaler, scale_zero,
                                                     store_dtype, store_nodata_value)

    undetect = scaler * (PD["data_undetect"] - scale_zero)

    metadata = {
        "nodata": store_nodata_value,
        "gain": 1./scaler,
        "offset": scale_zero,
        "undetect": undetect,
    }

    return prepared_forecast, metadata


def get_timesteps():
    """Return the nowcast timestep if it is regular"""
    runopt = PD["run_options"]
    ts = runopt.get("nowcast_timestep")
    if ts is None:
        return PD["data_source"]["timestep"]
    return ts


if __name__ == '__main__':
    run(test=True)
