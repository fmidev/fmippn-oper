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

    # generate suitable objects for passing to pysteps methods
    datasource, nowcast_kwargs = generate_pysteps_setup()

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
    print(run_options)

    observations, obs_metadata = read_observations(startdate, datasource, importer)

    motion_field = optflow(observations)

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
    else:
        ensemble_forecast = None
        ens_meta = dict()

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
    PD["ENSEMBLE_SIZE"] = PD.get("nowcast_options").get("n_ens_members")
    use_old_format = PD["output_options"].get("use_old_format", False)
    del PD["data_source"]
    del PD["data_options"]
    del PD["motion_options"]
    del PD["nowcast_options"]
    del PD["run_options"]
    del PD["output_options"]
    del PD["logging"]
    if store_meta["seed"] is None:  # Cannot write None to HDF5
        del store_meta["seed"]
    if "SEED" in PD:
        del PD["SEED"]

    # WRITE OUTPUT TO A FILE
    if use_old_format:
        write_to_file(startdate, gen_output, nc_fname, store_meta)
    else:
        # FIXME: Add new ODIM-style function here
        print("Placeholder: Write using new format")
    # FIXME: Disable temporarily, uncomment when write_to_file ignores dictionaries
    # ~ log("info", "Finished writing output to a file.")
    # ~ log("info", "Run complete. Exiting.")


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
        **kwargs -- additional keyword arguments passed to importer method

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
        **kwargs -- additional keyword arguments passed to optical flow method

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
        **kwargs -- additional keyword arguments passed to nowcast method

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
        **kwargs -- additional keyword arguments passed to nowcast method

    Output:
        function -- a function object

    Raise ValueError for invalid `module` selectors.
    """
    if module == "pysteps":
        return pysteps.nowcasts.get_method(PD["run_options"]["deterministic_method"], **kwargs)
    # Add more options here

    raise ValueError("Unknown module {} for deterministic method".format(module))

def generate_pysteps_setup():
    """Generate `datasource` and `nowcast_kwargs` objects that are suitable
    for using in pysteps nowcasting methods."""
    # Paths, importers etc.
    datasource = PD.get("data_source")
    # NOTE: This is for backwards compability, can be removed at some point
    if datasource is None:
        datasource = pystepsrc["data_sources"][PD["DOMAIN"]]
    datasource["root_path"] = os.path.expanduser(datasource["root_path"])

    # kwargs for nowcasting method
    nowcast_kwargs = PD.get("nowcast_options")

    # This threshold is used in masking and probability masking
    # rrate units need to be transformed to decibel, so that comparisons can be done
    r_thr = PD["RAIN_THRESHOLD"]
    log("debug", f"Using RAIN_THRESHOLD {r_thr} value as prob. match threshold")

    if PD["VALUE_DOMAIN"] == "rrate":
        nowcast_kwargs["R_thr"] = 10.0 * np.log10(r_thr)
    else:
        nowcast_kwargs["R_thr"] = r_thr

    return datasource, nowcast_kwargs

def read_observations(startdate, datasource, importer):
    """Read observations from archives using pysteps methods. Also threshold
    the input data and (optionally) convert dBZ -> dBR based on configuration
    parameters."""
    try:
        filelist = pysteps.io.find_by_date(startdate,
                                           datasource["root_path"],
                                           datasource["path_fmt"],
                                           datasource["fn_pattern"],
                                           datasource["fn_ext"],
                                           datasource["timestep"],
                                           num_prev_files=PD["NUM_PREV_OBSERVATIONS"])
    except OSError as pysteps_error:
        error_msg = "Failed to read input data!"
        log("error", f"OSError was raised: {error_msg}")
        # Re-raise so traceback is shown in stdout and program stops
        raise OSError(error_msg) from pysteps_error

    # PGM files contain dBZ values
    obs, _, metadata = pysteps.io.readers.read_timeseries(filelist,
                                                          importer,
                                                          **datasource["importer_kwargs"])

    if PD["VALUE_DOMAIN"] == "rrate":
        obs, metadata = dbz_to_rrate(obs, metadata)

    obs, metadata = thresholding(obs, metadata, threshold=PD["RAIN_THRESHOLD"],
                                 norain_value=PD["NORAIN_VALUE"])

    if PD["VALUE_DOMAIN"] == "rrate":
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

def thresholding(data, metadata, threshold=None, norain_value=None, fill_nan=True):
    if threshold is None:
        threshold = PD["RAIN_THRESHOLD"]
    if norain_value is None:
        norain_value = PD["NORAIN_VALUE"]

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

    if PD["FIELD_VALUES"] == "dbz" and metadata["unit"] == "mm/h":
        forecast, meta = rrate_to_dbz(forecast, meta)

    if PD["FIELD_VALUES"] == "rrate" and metadata["unit"] == "dBZ":
        forecast, meta = dbz_to_rrate(forecast, meta)

    log("debug", f"{metadata['unit']}")
    norain_value = -32 if metadata["unit"] == "dBZ" else 0
    forecast, meta = thresholding(forecast, meta, norain_value=norain_value, fill_nan=False)

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
    scaler = PD["SCALER"]
    scale_zero = PD["SCALE_ZERO"]
    if scale_zero in [None, "auto"]:
        scale_zero = np.nanmin(forecast)
    prepared_forecast = utils.prepare_fct_for_saving(forecast, scaler, scale_zero,
                                                     store_dtype, store_nodata_value)

    metadata = {
        "nodata": store_nodata_value,
        "gain": 1./scaler,
        "offset": scale_zero,
    }

    return prepared_forecast, metadata

def write_to_file(startdate, gen_output, nc_fname, metadata=None):
    """Write output to a HDF5 file.

    Input:
        startdate -- nowcast analysis time (datetime object)
        gen_output -- dictionary containing generated nowcasts
        nc_fname -- filename for output HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """
    ensemble_forecast = gen_output.get("ensemble_forecast", None)
    unperturbed = gen_output.get("unperturbed", None)
    deterministic = gen_output.get("deterministic", None)
    motion_field = gen_output.get("motion_field", None)
    ensemble_motion = gen_output.get("ensemble_motion", None)

    if metadata is None:
        metadata = dict()

    if all((dataset is None for dataset in gen_output.values())):
        print("Nothing to store")
        log("warning", "Nothing to store into .h5 file. Skipping.")
        return None

    ensemble_forecast, ens_scale_meta = prepare_data_for_writing(ensemble_forecast)
    deterministic, det_scale_meta = prepare_data_for_writing(deterministic)
    unperturbed, unpert_scale_meta = prepare_data_for_writing(unperturbed)

    with h5py.File(os.path.join(PD["OUTPUT_PATH"], nc_fname), 'w') as outf:
        if ensemble_forecast is not None and PD["STORE_ENSEMBLE"]:
            for eidx in range(PD["ENSEMBLE_SIZE"]):
                ens_grp = outf.create_group("member-{:0>2}".format(eidx))
                utils.store_timeseries(ens_grp,
                                       ensemble_forecast[eidx, :, :, :],
                                       startdate,
                                       timestep=PD["NOWCAST_TIMESTEP"],
                                       metadata=ens_scale_meta)

        if ensemble_motion is not None and PD["STORE_PERTURBED_MOTION"]:
            for eidx in range(PD["ENSEMBLE_SIZE"]):
                try:
                    ens_grp = outf["member-{:0>2}".format(eidx)]
                except KeyError:
                    ens_grp = outf.create_group("member-{:0>2}".format(eidx))
                ens_grp.create_dataset("motion", data=ensemble_motion[eidx])

        if unperturbed is not None and PD["STORE_UNPERTURBED"]:
            unpert_grp = outf.create_group("unperturbed")
            utils.store_timeseries(unpert_grp, unperturbed[0, :, :, :], startdate,
                                   timestep=PD["NOWCAST_TIMESTEP"],
                                   metadata=det_scale_meta)

        if deterministic is not None and PD["STORE_DETERMINISTIC"]:
            det_grp = outf.create_group("deterministic")
            utils.store_timeseries(det_grp, deterministic, startdate,
                                   timestep=PD["NOWCAST_TIMESTEP"],
                                   metadata=unpert_scale_meta)

        if PD["STORE_MOTION"]:
            outf.create_dataset("motion", data=motion_field)

        meta = outf.create_group("meta")
        meta.attrs["nowcast_started"] = dt.datetime.strftime(metadata["time_at_start"],
                                                             PD["OUTPUT_TIME_FORMAT"])
        meta.attrs["nowcast_ended"] = dt.datetime.strftime(metadata["time_at_end"],
                                                           PD["OUTPUT_TIME_FORMAT"])
        meta.attrs["nowcast_units"] = metadata.get("unit", "Unknown")
        meta.attrs["nowcast_seed"] = metadata.get("seed", "Unknown")
        meta.attrs["nowcast_init_time"] = dt.datetime.strftime(startdate, "%Y%m%d%H%M")

        pd_meta = meta.create_group("configuration")
        for key, value in PD.items():
            if key in ["LOG_FOLDER", "LOG_LEVEL", "OUTPUT_PATH"]:
                continue
            pd_meta.attrs[key] = str(value)

        proj_meta = meta.create_group("projection")
        for key, value in metadata["projection"].items():
            proj_meta.attrs[key] = value

    return None

if __name__ == '__main__':
    run(test=True)
