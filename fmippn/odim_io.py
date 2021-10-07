"""Writing functions for storing the PPN output in HDF5 files"""
import os

import h5py

import utils
from ppn_config import defaults

def write_deterministic_to_file(configuration, nowcast_data, filename=None, metadata=None):
    """Write deterministic output in ODIM HDF5 format..

    Input:
        configuration -- Object containing configuration parameters
        nowcast_data -- generated deterministic nowcast
        filename -- filename for output deterministic HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """
    #Output filename
    if filename is None:
        filename = os.path.join(defaults["output_options"]["path"], "deterministic.h5")

    _write(nowcast_data, filename, metadata, configuration=configuration, optype="det")


def write_ensemble_to_file(configuration, nowcast_data, filename=None, metadata=None):
    """Write ensemble output in ODIM HDF5 format..

    Input:
        configuration -- Object containing configuration parameters
        nowcast_data -- generated ensemble nowcasts
        filename -- filename for output ensemble HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """
    #Output filename
    if filename is None:
        filename = os.path.join(defaults["output_options"]["path"], "ensemble.h5")

    _write(nowcast_data, filename, metadata,configuration=configuration,  optype="ens")


def write_motion_to_file(configuration, motion_data, filename=None, metadata=None):
    """Write motion field output in ODIM HDF5 format..

    Input:
        configuration -- Object containing configuration parameters
        motion_data -- motion field data
        filename -- filename for output motion HDF5 file
        metadata -- dictionary containing nowcast metadata (optional)
    """
    #Output filename
    if filename is None:
        filename = os.path.join(defaults["output_options"]["path"], "motion.h5")

    _write(motion_data, filename, metadata,configuration=configuration,  optype="mot")

# FIXME: This logic should be converted to use a list of leadtimes instead of assuming regular timestep
def get_timesteps(configuration):
    """Return the nowcast timestep if it is regular"""
    runopt = configuration["run_options"]
    ts = runopt.get("nowcast_timestep")
    if ts is None:
        return configuration["data_source"]["timestep"]
    return ts

def _convert_motion_units(data_pxts, kmperpixel=1.0, timestep=1.0):
    """Convert atmospheric motion vectors from pixel/timestep units to m/s.

    Input:
        data_pxts -- motion vectors in "pixels per timestep" units
        kmperpixel -- kilometers in pixel
        timestep -- timestep lenght in minutes

    Output:
        data_ms -- motion vectors in m/s units
    """
    meters_per_pixel = kmperpixel * 1000
    seconds_in_timestep = timestep * 60
    # data unit conversion logic:
    # px_per_timestep = px_per_s * s_per_timestep
    #                 = px_per_km * km_per_s * s_per_timestep
    #                 = px_per_km * km_per_m * m_per_s * s_per_timestep
    # Solve for m_per_s:
    # m_per_s = px_per_timestep / (px_per_km * km_per_m * s_per_timestep)
    #         = px_per_timestep * km_per_px * m_per_km / s_per_timestep
    #         = px_per_timestep * meters_per_pixel / seconds_in_timestep
    data_ms = data_pxts * meters_per_pixel / seconds_in_timestep
    return data_ms

def _write(data, filename, metadata, configuration, optype=None):
    # Necessary input value checks, exit early if no need to store anything
    if data is None:
        print("Nothing to store")
        return None

    if optype not in {"det", "ens", "mot"}:
        print("Missing logic for this optype:", optype)
        return None

    if metadata is None:
        metadata = dict()
    else:
        metadata = metadata.copy()  # Do not modify original

    # Variables needed for nowcasts
    if optype in {"det", "ens"}:
        nowcast_timestep = get_timesteps(configuration)
        scale_meta = metadata.get("scale_meta", dict())

        startdate = metadata.get("startdate")
        if startdate is None:
            raise ValueError("missing startdate information from metadata dictionary")

    # Conversion of motion vector units from pixel/timestep -> m/s
    if optype == 'mot':
        motion_timestep = configuration["data_source"].get("timestep")
        motion_pixelsize = configuration.get("nowcast_options").get("kmperpixel")
        data = _convert_motion_units(
            data_pxts=data,
            kmperpixel=motion_pixelsize,
            timestep=motion_timestep
        )

    # Remove things that break HDF5 output from metadata dictionary
    if "scale_meta" in metadata:
        del metadata["scale_meta"]
    if "startdate" in metadata:
        del metadata["startdate"]
    if "seed" in metadata and metadata["seed"] is None:
        del metadata["seed"]

    with h5py.File(filename, 'w') as outf:
        # Initialize output file
        # Copy attribute groups /what, /where and /how from input to output
        utils.copy_odim_attributes(configuration["odim_metadata"], outf)

        # Common groups
        how_grp = outf["how"]
        how_grp.attrs["domain"] = configuration["nowcast_options"]["domain"]

        #Write AMVU and AMVV datasets and add attributes
        if optype == "mot":
            AMVU = data[0]
            AMVV = data[1]

            how_attrs = {
                "input_interval": 60 * motion_timestep,
                "kmperpixel": motion_pixelsize,
                "units": "m/s,"
            }

            amvu_grp = outf.create_group("/dataset1/data1")
            amvu_grp.create_dataset("data", data=AMVU)
            amvu_what_grp = amvu_grp.create_group("what")
            amvu_what_grp.attrs["quantity"] = "AMVU"
            amvu_how_grp = amvu_grp.create_group("how")
            for key, value in how_attrs.items():
                amvu_how_grp.attrs[key] = value

            amvv_grp = outf.create_group("/dataset1/data2")
            amvv_grp.create_dataset("data", data=AMVV)
            amvv_what_grp = amvv_grp.create_group("what")
            amvv_what_grp.attrs["quantity"] = "AMVV"
            amvv_how_grp = amvv_grp.create_group("how")
            for key, value in how_attrs.items():
                amvv_how_grp.attrs[key] = value

        #Write deterministic forecast timeseries in ODIM format
        elif optype == "det":
            for index in range(data.shape[0]):
                dset_grp=outf.create_group(f"/dataset{index+1}")

                #Add attributes to each dataset
                utils.store_odim_dset_attrs(dset_grp, index, startdate, nowcast_timestep)

                #Store data
                ts_point = data[index, :, :]
                data_grp=dset_grp.create_group("data1")
                data_grp.create_dataset("data",data=ts_point)

                #Store data/what group attributes
                utils.store_odim_data_what_attrs(data_grp, metadata, scale_meta)

            #Store PPN specific metadata into /how group
            how_grp.attrs["zr_a"] = configuration["data_options"]["zr_a"]
            how_grp.attrs["zr_b"] = configuration["data_options"]["zr_b"]
            # FIXME: "leadtimes" can be a list (irregular timesteps) -> take that into account
            how_grp.attrs["num_timesteps"] = configuration["run_options"]["leadtimes"]
            # FIXME: "nowcast_timestep" might not be constants, see above
            how_grp.attrs["nowcast_timestep"] = nowcast_timestep
            how_grp.attrs["max_leadtime"] = configuration["run_options"]["max_leadtime"]
            default_cascade_levels = defaults["nowcast_options"]["n_cascade_levels"]
            how_grp.attrs["n_cascade_levels"] = configuration["nowcast_options"].get("n_cascade_levels",
                                                                                     default_cascade_levels)

        #Write ensemble forecast timeseries in ODIM format
        elif optype == "ens":
            for index in range(data.shape[1]):
                dset_grp=outf.create_group(f"/dataset{index+1}")

                #Add attributes to each dataset
                utils.store_odim_dset_attrs(dset_grp, index, startdate, nowcast_timestep)

                # Store ensemble members
                for eidx in range(configuration["ensemble_size"]):

                    #Store data
                    ts_point = data[eidx, index, :, :]
                    data_grp=dset_grp.create_group(f"data{eidx+1}")
                    data_grp.create_dataset("data",data=ts_point)

                    #Store data/what group attributes
                    utils.store_odim_data_what_attrs(data_grp, metadata, scale_meta)

            #Store PPN specific metadata into /how group
            how_grp.attrs["zr_a"] = configuration["data_options"]["zr_a"]
            how_grp.attrs["zr_b"] = configuration["data_options"]["zr_b"]
            how_grp.attrs["seed"] = metadata.get("seed","Unknown")
            how_grp.attrs["ensemble_size"] = configuration["ensemble_size"]
            # FIXME: "leadtimes" can be a list (irregular timesteps) -> take that into account
            how_grp.attrs["num_timesteps"] = configuration["run_options"]["leadtimes"]
            # FIXME: "nowcast_timestep" might not be constants, see above
            how_grp.attrs["nowcast_timestep"] = nowcast_timestep
            how_grp.attrs["max_leadtime"] = configuration["run_options"]["max_leadtime"]
            default_cascade_levels = defaults["nowcast_options"]["n_cascade_levels"]
            how_grp.attrs["n_cascade_levels"] = configuration["nowcast_options"].get("n_cascade_levels",
                                                                                     default_cascade_levels)

    return None
