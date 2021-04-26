"""Main program for FMI-PPN.

FMI-PPN (Finnish Meteorological Institute Probabilistic Precipitation
Nowcaster) is a script for weather radar-based nowcasting. It is
heavily based on pySTEPS initiative.

For more information about pySTEPS, see https://pysteps.github.io .

Author: Petteri Karsisto
Year: 2019
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

    initialise_logging(log_folder=PD["logging"]["log_folder"],
                       log_fname="ppn-{:%Y%m%d}.log".format(dt.datetime.utcnow()))

    log("info", "Program starting.")
    log("debug", "Start setup")
    log("debug", f"Input: timestamp={timestamp}, config={config}")

    # Default is generate nowcasts from previous radar composite
    # We need to replace utcnow() minutes with nearest (floor) multiple of 5
    # to get a valid timestamp for the data
    # However, if a timestamp is given, that should be used.
    if timestamp is not None:
        startdate = dt.datetime.strptime(timestamp, "%Y%m%d%H%M")
    else:
        startdate = utils.utcnow_floored(increment=5)

    if nc_fname is None:
        nc_fname = "nc_{:%Y%m%d%H%M}.h5".format(startdate)

    # GENERAL SETUP

    # Paths, importers etc.
    datasource = PD.get("data_source")
    # NOTE: This is for backwards compability, can be removed at some point
    if datasource is None:
        datasource = pystepsrc["data_sources"][PD["DOMAIN"]]

    # Used methods
    importer = importer_method(name=datasource["importer"])
    optflow = optflow_method("pysteps")
    nowcaster = nowcast_method("pysteps")
    deterministic_nowcaster = deterministic_method("pysteps")

    log("debug", "Setup finished")

    # NOWCASTING

    log("info", "Generating nowcasts, starting from %s" % (startdate))

    time_at_start = dt.datetime.today()

    run_options = PD["run_options"]

    input_files = get_filelist(startdate, datasource)

    if datasource["importer"] in ["opera_hdf5", "odim_hdf5"]:
        input_quantity = datasource["importer_kwargs"]["qty"]
        odim_metadata = utils.get_odim_attrs_from_input(input_files[0][-1])  # input_files is a tuple of two lists
        data_undetect = utils.get_odim_data_undetect(input_files[0][-1],
                                                     datasource['importer_kwargs']['qty'])
    else:
        # Cannot read ODIM metadata from non-ODIM files (.pgm)
        odim_metadata = None
        data_undetect = -32
        input_quantity = "DBZH"
    PD["odim_metadata"] = odim_metadata
    PD["data_undetect"] = data_undetect
    PD["input_quantity"] = input_quantity

    # generate suitable objects for passing to pysteps methods
    if PD["run_options"]["nowcast_method"] in ["steps"]:
        nowcast_kwargs = generate_pysteps_setup()

    observations, obs_metadata = read_observations(input_files, datasource, importer)

    motion_field = optflow(observations, **PD.get("motion_options", dict()))

    # If seed is none, make a random seed.
    if PD["nowcast_options"].get("seed") is None:
        PD["nowcast_options"]["seed"] = random.randrange(2**32-1)

    # Regenerate ensemble motion
    if run_options.get("regenerate_perturbed_motion"):
        if PD["nowcast_options"].get("seed") is None:
            raise ValueError("Cannot regenerate motion field with unknown seed value!")
        log("info", "Regenerating ensemble motion fields...")
        ensemble_motion = regenerate_ensemble_motion(motion_field, nowcast_kwargs)
        log("info", "Finished regeneration.")
    else:
        ensemble_motion = None

    if run_options.get("run_ensemble"):
        ensemble_forecast, ens_meta = generate(observations, motion_field, nowcaster,
                                               nowcast_kwargs, metadata=obs_metadata)
        PD["ensemble_size"] = ensemble_forecast.shape[0]
    else:
        ensemble_forecast = None
        ens_meta = dict()
        PD["ensemble_size"] = None

    if run_options.get("run_deterministic"):
        deterministic, det_meta = generate_deterministic(observations[-1],
                                                         motion_field,
                                                         deterministic_nowcaster,
                                                         metadata=obs_metadata)
    else:
        deterministic = None
        det_meta = dict()

    time_at_end = dt.datetime.today()
    log("debug", "Finished nowcasting at %s" % time_at_end)
    log("info", "Finished nowcasting. Time elapsed: %s" % (time_at_end - time_at_start))

    gen_output = {
        "motion_field": motion_field,
        "ensemble_motion": ensemble_motion,
        "ensemble_forecast": ensemble_forecast,
        "deterministic": deterministic,
    }

    # Metadata for storage
    if "unit" in ens_meta:
        unit = ens_meta["unit"]
    elif "unit" in det_meta:
        unit = det_meta["unit"]
    else:
        unit = "Unknown"
    store_meta = {
        "unit": unit,
        "seed": PD["nowcast_options"]["seed"],
        "projection": {
            "projstr": obs_metadata["projection"],
            "x1": obs_metadata["x1"],
            "x2": obs_metadata["x2"],
            "y1": obs_metadata["y1"],
            "y2": obs_metadata["y2"],
            "xpixelsize": obs_metadata["xpixelsize"],
            "ypixelsize": obs_metadata["ypixelsize"],
            "origin": "upper",
        },
        "time_at_start": time_at_start,
        "time_at_end": time_at_end,
    }

    # FIXME: temporary hack to prevent crashes during saving the output
    # Remove these when the write_to_file function has been rewritten
    if store_meta["seed"] is None:  # Cannot write None to HDF5
        del store_meta["seed"]

    # WRITE OUTPUT TO A FILE
    if PD["output_options"].get("use_old_format", False):
        write_to_file(startdate, gen_output, nc_fname, store_meta)
    else:
        # Test writing ODIM output
        nc_det_fname=None
        nc_ens_fname=None
        nc_mot_fname=None
        write_odim_deterministic_to_file(startdate, gen_output, nc_det_fname, store_meta)
        write_odim_ensemble_to_file(startdate, gen_output, nc_ens_fname, store_meta)
        write_odim_motion_to_file(startdate, gen_output, nc_mot_fname, store_meta)

    log("info", "Finished writing output to a file.")
    log("info", "Run complete. Exiting.")


def initialise_logging(log_folder='./', log_fname='ppn.log'):
    """Wrapper for ppn_logger.config_logging() method. Does nothing if writing
    to log is not enabled."""
    if PD["logging"]["write_log"]:
        full_path = os.path.expanduser(log_folder)
        ppn_logger.config_logging(os.path.join(full_path, log_fname),
                                  level=PD["logging"]["log_level"])

def log(level, msg, *args, **kwargs):
    """Wrapper for ppn_logger. Function does nothing if writing to log is
    not enabled."""
    if PD["logging"]["write_log"]:
        ppn_logger.write_to_log(level, msg, *args, **kwargs)

def importer_method(module="pysteps", **kwargs):
    """Wrapper for easily switching between modules which provide data importer
    methods.

    Input:
        module -- parameter for if/else block (default="pysteps")
        **kwargs -- additional keyword arguments passed to importer method getter

    Output:
        function -- a function object

    Raise ValueError for invalid `module` selectors.
    """
    if module == "pysteps":
        return pysteps.io.get_method(method_type="importer", **kwargs)
    # Add more options here

    raise ValueError("Unknown module {} for importer method".format(module))

def optflow_method(module="pysteps", **kwargs):
    """Wrapper for easily switching between modules which provide optical flow
    methods.

    Input:
        module -- parameter for if/else block (default="pysteps")
        **kwargs -- additional keyword arguments passed to optical flow method getter

    Output:
        function -- a function object

    Raise ValueError for invalid `module` selectors.
    """
    if module == "pysteps":
        return pysteps.motion.get_method(PD["run_options"]["motion_method"], **kwargs)
    # Add more options here

    raise ValueError("Unknown module {} for optical flow method".format(module))

def nowcast_method(module="pysteps", **kwargs):
    """Wrapper for easily switching between modules which provide nowcasting
    methods.

    Input:
        module -- parameter for if/else block (default="pysteps")
        **kwargs -- additional keyword arguments passed to nowcast method getter

    Output:
        function -- a function object

    Raise ValueError for invalid `module` selectors.
    """
    if module == "pysteps":
        return pysteps.nowcasts.get_method(PD["run_options"]["nowcast_method"], **kwargs)
    # Add more options here

    raise ValueError("Unknown module {} for nowcast method".format(module))

def deterministic_method(module="pysteps", **kwargs):
    """Wrapper for easily switching between modules which provide deterministic
    nowcasting methods.

    Input:
        module -- parameter for if/else block (default="pysteps")
        **kwargs -- additional keyword arguments passed to nowcast method getter

    Output:
        function -- a function object

    Raise ValueError for invalid `module` selectors.
    """
    if module == "pysteps":
        return pysteps.nowcasts.get_method(PD["run_options"]["deterministic_method"], **kwargs)
    # Add more options here

    raise ValueError("Unknown module {} for deterministic method".format(module))

def generate_pysteps_setup():
    """Generate `nowcast_kwargs` objects that are suitable
    for using in pysteps nowcasting methods."""
    # kwargs for nowcasting method
    nowcast_kwargs = PD.get("nowcast_options")

    # This threshold is used in masking and probability masking
    # rrate units need to be transformed to decibel, so that comparisons can be done
    # Check if forecast is done for different quantity than input and convert if necessary
    r_thr = PD["data_options"].get("rain_threshold")

    fct_qty = PD["run_options"].get("forecast_as_quantity", PD["input_quantity"])
    log("debug", f"Using rain_threshold={r_thr} as prob. match threshold")

    zr_a = PD["data_options"]["zr_a"]
    zr_b = PD["data_options"]["zr_b"]

    if PD["input_quantity"] in ["DBZH"] and fct_qty in ["rrate", "RATE"]:
        r_thr = (r_thr / zr_a) ** (1. / zr_b)
        nowcast_kwargs["R_thr"] = max(10.0 * np.log10(r_thr), 0)  #
        log("info", 'Converted RATE rain_threshold to decibel units ("dBR").')
    elif PD["input_quantity"] in ["RATE"] and fct_qty in ["dbz", "DBZH"]:
        r_thr = zr_a * r_thr ** zr_b
        nowcast_kwargs["R_thr"] = r_thr
    else:
        nowcast_kwargs["R_thr"] = r_thr

    PD["converted_rain_thr"] = r_thr  # DBZH or non-decibel RATE is used in thresholding

    return nowcast_kwargs

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

def read_observations(filelist, datasource, importer):
    """Read observations from archives using pysteps methods. Also threshold
    the input data and (optionally) convert dBZ -> dBR based on configuration
    parameters."""
    # PGM files contain dBZ values
    obs, _, metadata = pysteps.io.readers.read_timeseries(filelist,
                                                          importer,
                                                          **datasource["importer_kwargs"])

    fct_qty = PD["run_options"].get("forecast_as_quantity", PD["input_quantity"])
    if PD["input_quantity"] in ["DBZH"] and fct_qty in ["rrate", "RATE"]:
        obs, metadata = dbz_to_rrate(obs, metadata)
    elif PD["input_quantity"] in ["RATE"] and fct_qty in ["dbz", "DBZH"]:
        obs, metadata = rrate_to_dbz(obs, metadata)

    obs, metadata = thresholding(obs, metadata, threshold=PD["converted_rain_thr"],
                                 norain_value=PD["run_options"]["steps_set_no_rain_to_value"])

    if fct_qty in ["RATE", "rrate"]:
        obs, metadata = transform_to_decibels(obs, metadata)

    return obs, metadata

def dbz_to_rrate(data, metadata):
    return pysteps.utils.conversion.to_rainrate(data, metadata, PD["data_options"]["zr_a"],
                                                PD["data_options"]["zr_b"])

def rrate_to_dbz(data, metadata):
    return pysteps.utils.conversion.to_reflectivity(data, metadata, PD["data_options"]["zr_a"],
                                                    PD["data_options"]["zr_b"])

def transform_to_decibels(data, metadata, inverse=False):
    """Transform data to decibel units. Assumes thresholded data.

    If argument `inverse` is True, perform transform from decibel units to
    normal units.
    """
    threshold = metadata["threshold"]
    zerovalue = metadata["zerovalue"]
    transform = metadata["transform"]

    data = data.copy()
    metadata = metadata.copy()

    if inverse:
        log("debug", "Converting data FROM dB units")
        new_transform = None
        new_threshold = 10.0 ** (threshold / 10.0)
        new_zerovalue = 10.0 ** (zerovalue / 10.0)
        data = 10.0 ** (data / 10.0)
    else:
        log("debug", "Converting data TO dB units")
        new_transform = "dB"
        new_threshold = 10.0 * np.log10(threshold)
        new_zerovalue = 10.0 * np.log10(zerovalue)
        data = 10.0 * np.log10(data)

    metadata["transform"] = new_transform
    metadata["zerovalue"] = new_zerovalue
    metadata["threshold"] = new_threshold

    return data, metadata

def thresholding(data, metadata, threshold, norain_value, fill_nan=True):
    """Set values under 'threshold' to 'norain_value'. Optionally replace np.nan with 'norain_value'.
    """
    if fill_nan:
        data[~np.isfinite(data)] = norain_value
    data[data < threshold] = norain_value

    metadata["zerovalue"] = norain_value
    metadata["threshold"] = threshold

    return data, metadata

def generate(observations, motion_field, nowcaster, nowcast_kwargs, metadata=None):
    """Generate ensemble nowcast using pysteps nowcaster."""
    forecast = nowcaster(observations, motion_field, PD["run_options"]["leadtimes"],
                         **nowcast_kwargs)

    if (metadata["unit"] == "mm/h") and (metadata["transform"] == "dB"):
        forecast, meta = transform_to_decibels(forecast, metadata, inverse=True)
    else:
        meta = metadata

    # FIXME: Logic is probably unnecessarily convoluted, needs simplifying and probably reordering
    # TODO: Should probably move this to own function...
    # Quantity conversion from calculation quantity to output quantity, if they are different
    out_qty = PD["output_options"].get("as_quantity", None)
    if out_qty is None:
        out_qty = PD["input_quantity"]

    if out_qty in ["dbz", "DBZH"] and metadata["unit"] == "mm/h":
        forecast, meta = rrate_to_dbz(forecast, meta)

    elif out_qty in ["rrate", "RATE"] and metadata["unit"] == "dBZ":
        forecast, meta = dbz_to_rrate(forecast, meta)

    # Set values under rain threshold to original undetect value
    # But first, check that the units are correct and convert data_undetect to other units if needed
    rain_threshold = PD["data_options"].get("rain_threshold")
    norain_for_output = PD["data_undetect"]
    if PD["input_quantity"] in ["dbz", "DBZH"] and out_qty in ["rrate", "RATE"]:
        # R = (Z / zr_a) ** (1.0 / zr_b)
        zr_a = PD["data_options"]["zr_a"]
        zr_b = PD["data_options"]["zr_b"]
        norain_for_output = (norain_for_output / zr_a) ** (1. / zr_b)
        rain_threshold = (rain_threshold / zr_a) ** (1. / zr_b)
    elif PD["input_quantity"] in ["rrate", "RATE"] and out_qty in ["dbz", "DBZH"]:
        # Z = zr_a * R ** zr_b
        zr_a = PD["data_options"]["zr_a"]
        zr_b = PD["data_options"]["zr_b"]
        norain_for_output = zr_a * norain_for_output ** zr_b
        rain_threshold = zr_a * rain_threshold ** zr_b

    forecast, meta = thresholding(forecast, meta, threshold=rain_threshold,
                                  norain_value=norain_for_output, fill_nan=False)
    PD["out_rain_threshold"] = rain_threshold
    PD["out_norain_value"] = norain_for_output
    if meta is None:
        meta = dict()
    return forecast, meta

def generate_deterministic(observations, motion_field, nowcaster, nowcast_kwargs=None,
                           metadata=None):
    """Generate a deterministic nowcast using semilagrangian extrapolation"""
    # Extrapolation scheme doesn't use the same nowcast_kwargs as steps
    if nowcast_kwargs is None:
        nowcast_kwargs = dict()
    forecast, meta = generate(observations, motion_field, nowcaster, nowcast_kwargs,
                              metadata)
    return forecast, meta

def regenerate_ensemble_motion(motion_field, nowcast_kwargs):
    """Generate motion perturbations the same way as pysteps.nowcasts.steps function.

    This is a workaround for obtaining perturbed motion fields from steps
    calculations, as pysteps doesn't currently give them as output. (2019-09-02)
    """
    pixelsperkm = 1./nowcast_kwargs["kmperpixel"]
    timestep = nowcast_kwargs["timestep"]
    pert_params = nowcast_kwargs["vel_pert_kwargs"]

    # (edited from pysteps.nowcasts.steps.forecast function)
    # initialize the random generators
    seed = nowcast_kwargs["seed"]
    if nowcast_kwargs["vel_pert_method"] is not None:
        randgen_prec = []
        randgen_motion = []
        np.random.seed(seed)
        for _ in range(nowcast_kwargs["n_ens_members"]):
            new_state = np.random.RandomState(seed)
            randgen_prec.append(new_state)
            seed = new_state.randint(0, high=1e9)
            new_state = np.random.RandomState(seed)
            randgen_motion.append(new_state)
            seed = new_state.randint(0, high=1e9)
    # (copypaste ends here)

    ensemble_motions = []
    for i, random_state in enumerate(randgen_motion):
        init_perturbations = pysteps.noise.motion.initialize_bps(motion_field,
                                                                 pixelsperkm,
                                                                 timestep,
                                                                 p_par=pert_params["p_par"],
                                                                 p_perp=pert_params["p_perp"],
                                                                 randstate=random_state)
        perturbations = pysteps.noise.motion.generate_bps(init_perturbations, timestep*(i+1))
        perturbed = motion_field + perturbations
        ensemble_motions.append(perturbed)

    return ensemble_motions

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

# FIXME: This logic should be converted to use a list of leadtimes instead of assuming regular timestep
def get_timesteps():
    """Return the nowcast timestep if it is regular"""
    runopt = PD["run_options"]
    ts = runopt.get("nowcast_timestep")
    if ts is None:
        return PD["data_source"]["timestep"]
    return ts

def write_to_file(startdate, gen_output, nc_fname, metadata=None):
    """Write output to a HDF5 file.

    Input:
        startdate -- nowcast analysis time (datetime object)
        gen_output -- dictionary containing generated nowcasts
        nc_fname -- filename for output HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """
    ensemble_forecast = gen_output.get("ensemble_forecast", None)
    deterministic = gen_output.get("deterministic", None)
    motion_field = gen_output.get("motion_field", None)
    ensemble_motion = gen_output.get("ensemble_motion", None)

    if metadata is None:
        metadata = dict()

    if all((dataset is None for dataset in gen_output.values())):
        print("Nothing to store")
        log("warning", "Nothing to store into .h5 file. Skipping.")
        return None

    output_options = PD["output_options"]
    # We don't support irregular nowcast outputs here
    nowcast_timestep = get_timesteps()

    ensemble_forecast, ens_scale_meta = prepare_data_for_writing(ensemble_forecast)
    deterministic, det_scale_meta = prepare_data_for_writing(deterministic)

    with h5py.File(os.path.join(output_options["path"], nc_fname), 'w') as outf:
        if ensemble_forecast is not None and output_options["store_ensemble"]:
            for eidx in range(PD["ensemble_size"]):
                ens_grp = outf.create_group("member-{:0>2}".format(eidx))
                utils.store_timeseries(ens_grp,
                                       ensemble_forecast[eidx, :, :, :],
                                       startdate,
                                       timestep=nowcast_timestep,
                                       metadata=ens_scale_meta)

        if ensemble_motion is not None and output_options["store_perturbed_motion"]:
            for eidx in range(PD["ensemble_size"]):
                try:
                    ens_grp = outf["member-{:0>2}".format(eidx)]
                except KeyError:
                    ens_grp = outf.create_group("member-{:0>2}".format(eidx))
                ens_grp.create_dataset("motion", data=ensemble_motion[eidx])

        if deterministic is not None and output_options["store_deterministic"]:
            det_grp = outf.create_group("deterministic")
            utils.store_timeseries(det_grp, deterministic, startdate,
                                   timestep=nowcast_timestep,
                                   metadata=det_scale_meta)

        if output_options["store_motion"]:
            outf.create_dataset("motion", data=motion_field)

        meta = outf.create_group("meta")
        # configuration "OUTPUT_TIME_FORMAT" is removed, new output uses ODIM standard
        meta.attrs["nowcast_started"] = dt.datetime.strftime(metadata["time_at_start"],
                                                             "%Y-%m-%d %H:%M:%S")
        meta.attrs["nowcast_ended"] = dt.datetime.strftime(metadata["time_at_end"],
                                                           "%Y-%m-%d %H:%M:%S")
        meta.attrs["nowcast_units"] = metadata.get("unit", "Unknown")
        meta.attrs["nowcast_seed"] = metadata.get("seed", "Unknown")
        meta.attrs["nowcast_init_time"] = dt.datetime.strftime(startdate, "%Y%m%d%H%M")

        # Old configurations - may be used by postprocessing scripts
        old_style_configs = {
            # Method selections
            #"DOMAIN": "fmi", # postprocessing defines this instead of reading it here
            "VALUE_DOMAIN": "rrate" if PD["run_options"]["forecast_as_quantity"] == "RATE" else "dbz",  # Unused?
            # Z-R conversion parameters
            "ZR_A": PD["data_options"]["zr_a"],  #
            "ZR_B": PD["data_options"]["zr_b"],  #
            # Nowcasting parameters
            "NOWCAST_TIMESTEP": nowcast_timestep,  #
            "MAX_LEADTIME": PD["run_options"]["max_leadtime"],  #
            "NUM_TIMESTEPS": PD["run_options"]["leadtimes"],  #
            "ENSEMBLE_SIZE": PD["ensemble_size"],  #
            "NUM_CASCADES": PD["nowcast_options"].get("n_cascade_levels", 6),  # Unused?
            "RAIN_THRESHOLD": PD["out_rain_threshold"],  # Unused?
            "NORAIN_VALUE": PD["out_norain_value"],  #
            "KMPERPIXEL": PD["nowcast_options"]["kmperpixel"],  # Unused?
            "CALCULATION_DOMAIN": PD["nowcast_options"]["domain"],  # Unused?
            "VEL_PERT_KWARGS": PD["nowcast_options"]["vel_pert_kwargs"],
            # Storing parameters
            "FIELD_VALUES": PD["output_options"]["as_quantity"],  # Unused?
            "STORE_DETERMINISTIC": output_options["store_deterministic"],  #
            "STORE_PERTURBED_MOTION": output_options["store_perturbed_motion"],  #
        }

        pd_meta = meta.create_group("configuration")
        for key, value in old_style_configs.items():
            pd_meta.attrs[key] = str(value)

        proj_meta = meta.create_group("projection")
        for key, value in metadata["projection"].items():
            proj_meta.attrs[key] = value

    return None


def write_odim_deterministic_to_file(startdate, gen_output, nc_det_fname=None, metadata=None):
    """Write deterministic output in ODIM HDF5 format..

    Input:
        startdate -- nowcast analysis time (datetime object)
        gen_output -- dictionary containing generated nowcasts
        nc_det_fname -- filename for output deterministic HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """

    if metadata is None:
        metadata = dict()

    deterministic = gen_output.get("deterministic", None)

    #Output filename
    if nc_det_fname is None:
        nc_det_fname = "nc_det_{:%Y%m%d%H%M}.h5".format(startdate)

    if all((dataset is None for dataset in gen_output.values())):
        print("Nothing to store")
        log("warning", "Nothing to store into .h5 file. Skipping.")
        return None

    deterministic, det_scale_meta = prepare_data_for_writing(deterministic)

    #Write deterministic forecast in ODIM format
    if deterministic is not None and PD["output_options"]["store_deterministic"]:
        with h5py.File(os.path.join(PD["output_options"]["path"], nc_det_fname), 'w') as outf:

            #Copy attribute groups /what, /where and /how from input to output
            #utils.copy_odim_attributes(infile,outf)
            utils.copy_odim_attributes(PD["odim_metadata"],outf)

            # Write timeseries
            nowcast_timestep = get_timesteps()
            for index in range(deterministic.shape[0]):
                timestep=nowcast_timestep
                dset_grp=outf.create_group(f"/dataset{index+1}")

                #Add attributes to each dataset
                utils.store_odim_dset_attrs(dset_grp, index, startdate, timestep)

                #Store data
                ts_point = deterministic[index, :, :]
                data_grp=dset_grp.create_group("data1")
                data_grp.create_dataset("data",data=ts_point)

                #Store data/what group attributes
                utils.store_odim_data_what_attrs(data_grp,metadata,det_scale_meta)

                #Store PPN specific metadata into /how group
                how_grp=outf["how"]
                how_grp.attrs["zr_a"] = PD["data_options"]["zr_a"]
                how_grp.attrs["zr_b"] = PD["data_options"]["zr_b"]
                how_grp.attrs["domain"] = PD["nowcast_options"]["domain"]
                how_grp.attrs["num_timesteps"] = PD["run_options"]["leadtimes"]
                how_grp.attrs["nowcast_timestep"] = timestep
                how_grp.attrs["n_cascade_levels"] = PD["nowcast_options"].get("n_cascade_levels", 6)
                how_grp.attrs["max_leadtime"] = PD["run_options"]["max_leadtime"]


    return None



def write_odim_motion_to_file(startdate, gen_output, nc_mot_fname=None, metadata=None):
    """Write motion field output in ODIM HDF5 format..

    Input:
        startdate -- nowcast analysis time (datetime object)
        gen_output -- dictionary containing generated nowcasts
        nc_mot_fname -- filename for output motion HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """

    if metadata is None:
        metadata = dict()

    motion_field = gen_output.get("motion_field", None)

    #Output filename
    if nc_mot_fname is None:
        nc_mot_fname = "nc_motion_{:%Y%m%d%H%M}.h5".format(startdate)

    if all((dataset is None for dataset in gen_output.values())):
        print("Nothing to store")
        log("warning", "Nothing to store into .h5 file. Skipping.")
        return None

    #Write motion field in ODIM format
    if PD["output_options"]["store_motion"]:
        with h5py.File(os.path.join(PD["output_options"]["path"], nc_mot_fname), 'w') as outf:

            #TBD! Change from pix/s to m/s
            AMVU=motion_field[0]
            AMVV=motion_field[1]

            #Write AMVU and AMVV datasets and add attributes
            amvu_grp=outf.create_group("/dataset1/data1")
            amvu_grp.create_dataset("data", data=AMVU)
            amvu_what_grp=amvu_grp.create_group("what")
            amvu_what_grp.attrs["quantity"]="AMVU"

            amvv_grp=outf.create_group("/dataset1/data2")
            amvv_grp.create_dataset("data", data=AMVV)
            amvv_what_grp=amvv_grp.create_group("what")
            amvv_what_grp.attrs["quantity"]="AMVV"

            #Copy attribute groups /what, /where and /how from input to output
            #utils.copy_odim_attributes(infile,outf)
            utils.copy_odim_attributes(PD["odim_metadata"],outf)

            #Store PPN specific metadata into /how group
            how_grp=outf["how"]
            how_grp.attrs["domain"] = PD["nowcast_options"]["domain"]
            how_grp.attrs["n_cascade_levels"] = PD["nowcast_options"].get("n_cascade_levels", 6)

    return None



def write_odim_ensemble_to_file(startdate, gen_output, nc_ens_fname=None, metadata=None):
    """Write ensemble output in ODIM HDF5 format..

    Input:
        startdate -- nowcast analysis time (datetime object)
        gen_output -- dictionary containing generated nowcasts
        nc_ens_fname -- filename for output ensemble HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """

    if metadata is None:
        metadata = dict()

    ensemble_forecast = gen_output.get("ensemble_forecast", None)

    #Output filename
    if nc_ens_fname is None:
        nc_ens_fname = "nc_ens_{:%Y%m%d%H%M}.h5".format(startdate)

    if all((dataset is None for dataset in gen_output.values())):
        print("Nothing to store")
        log("warning", "Nothing to store into .h5 file. Skipping.")
        return None

    ensemble_forecast, ens_scale_meta = prepare_data_for_writing(ensemble_forecast)

    #Write ensemble forecast in ODIM format
    if ensemble_forecast is not None and PD["output_options"]["store_ensemble"]:
        with h5py.File(os.path.join(PD["output_options"]["path"], nc_ens_fname), 'w') as outf:

            #Copy attribute groups /what, /where and /how from input to output
            #utils.copy_odim_attributes(infile,outf)
            utils.copy_odim_attributes(PD["odim_metadata"],outf)

            # Write timeseries
            nowcast_timestep = get_timesteps()
            for index in range(ensemble_forecast.shape[1]):
                timestep=nowcast_timestep
                dset_grp=outf.create_group(f"/dataset{index+1}")

                #Add attributes to each dataset
                utils.store_odim_dset_attrs(dset_grp, index, startdate, timestep)

                # Store ensemble members
                for eidx in range(PD["ensemble_size"]):

                    #Store data
                    ts_point = ensemble_forecast[eidx, index, :, :]
                    data_grp=dset_grp.create_group(f"data{eidx+1}")
                    data_grp.create_dataset("data",data=ts_point)

                    #Store data/what group attributes
                    utils.store_odim_data_what_attrs(data_grp,metadata,ens_scale_meta)

                    #Store PPN specific metadata into /how group
                    how_grp=outf["how"]
                    how_grp.attrs["zr_a"] = PD["data_options"]["zr_a"]
                    how_grp.attrs["zr_b"] = PD["data_options"]["zr_b"]
                    how_grp.attrs["seed"] = metadata.get("seed","Unknown")
                    how_grp.attrs["domain"] = PD["nowcast_options"]["domain"]
                    how_grp.attrs["ensemble_size"] = PD["ensemble_size"]
                    how_grp.attrs["num_timesteps"] = PD["run_options"]["leadtimes"]
                    how_grp.attrs["nowcast_timestep"] = timestep
                    how_grp.attrs["n_cascade_levels"] = PD["nowcast_options"].get("n_cascade_levels", 6)
                    how_grp.attrs["max_leadtime"] = PD["run_options"]["max_leadtime"]

    return None





if __name__ == '__main__':
    run(test=True)
