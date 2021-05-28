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

from eccodes import *

import matplotlib.pyplot as plt


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

    print('input_files:', input_files)
        
    # Read observations from grib files
    observations=[]
    for grib_file in input_files[0]:
        print("grib_file:",grib_file)
        temps, temps_min, temps_max, dtime, mask_nodata, nodata, longitudes, latitudes=read_grib(grib_file)
        temps = temps[0,:,:]
        print('temps_min_max:', temps_min, temps_max)
        observations.append(temps)

    observations = np.concatenate([temps[None, :, :] for temps in observations])
    
    #For vet only 2 or 3 observation fields can be used instead of 4 used by lucaskanade!
    if PD["run_options"]["motion_method"] == "vet":
        print("obs shape orig",observations.shape)
        observations = observations[2:4]
        print("obs shape new",observations.shape)
        
    print('observations:',observations)
    print('min observations:',np.min(observations))
    print('max observations',np.max(observations))
    
    motion_field = optflow(observations, **PD.get("motion_options", dict()))

    print('motion_field:', motion_field)
    print('motion_field min max', np.min(motion_field), np.max(motion_field))

    #Read POT field
    last_obs_file = input_files[0][3]
    obs_path = os.path.split(last_obs_file)[0]
    obs_file = os.path.split(last_obs_file)[1]
    pot_file = obs_path + '/pot_' + obs_file
    print('pot_file',pot_file)

    temps, temps_min, temps_max, dtime, mask_nodata, nodata, longitudes, latitudes=read_grib(pot_file)
    pot_observations = temps[0,:,:]

    #deterministic_forecast = generate_deterministic(observations[-1], motion_field, deterministic_nowcaster)
    deterministic_forecast = generate_deterministic(pot_observations, motion_field, deterministic_nowcaster)
    
    print('deterministic_forecast:',deterministic_forecast)
    print('deterministic_forecast.shape:',deterministic_forecast.shape)

#    for timestep in range(deterministic_forecast.shape[0]):
#        print(timestep)
#        fname = "testfig_timestep=" + str(timestep) + ".png"
#        data=deterministic_forecast[timestep,:,:]
#        plt.imshow(data,origin='upper')
#        plt.savefig(fname)
#        plt.show()

    #Write deterministic forecast to grib file
    if nc_fname is None:
        nc_fname = "pot_nowcast_{:%Y%m%d%H%M}.grib2".format(startdate)
        
    write_grib(startdate,deterministic_forecast,pot_file,nc_fname)
        
        
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


def read_grib(image_grib_file,added_hours=0):

    dtime = []
    tempsl = []
    latitudes = []
    longitudes = []

    with GribFile(image_grib_file) as grib:
        for msg in grib:
            ni = msg["Ni"]
            nj = msg["Nj"]
            forecast_time = dt.datetime.strptime("{:d}/{:02d}".format(msg["dataDate"], int(msg["dataTime"]/100)), "%Y%m%d/%H") + dt.timedelta(hours=msg["forecastTime"])
            dtime.append(forecast_time)
            tempsl.append(np.asarray(msg["values"]).reshape(nj, ni))
            latitudes.append(np.asarray(msg["latitudes"]).reshape(nj, ni))
            longitudes.append(np.asarray(msg["longitudes"]).reshape(nj, ni))
    temps = np.asarray(tempsl)
    latitudes = np.asarray(latitudes)
    longitudes = np.asarray(longitudes)
    latitudes = latitudes[0,:,:]
    longitudes = longitudes[0,:,:]
    nodata = 9999
    mask_nodata = np.ma.masked_where(temps == nodata,temps)
    if len(temps[np.where(~np.ma.getmask(mask_nodata))])>0:
        temps_min = temps[np.where(~np.ma.getmask(mask_nodata))].min()
        temps_max = temps[np.where(~np.ma.getmask(mask_nodata))].max()
    else:
        print("input " + image_grib_file + " contains only missing data!")
        temps_min = nodata
        temps_max = nodata
    if type(dtime) == list:
        dtime = [(i+dt.timedelta(hours=added_hours)) for i in dtime]
    else:
        dtime = dtime+dt.timedelta(hours=added_hours)
    return temps, temps_min, temps_max, dtime, mask_nodata, nodata, longitudes, latitudes



def write_grib(startdate,deterministic_forecast,image_grib_file,nc_fname):
    ''' (Almost) all the metadata is copied from modeldata.grib2
     deterministic_forecast: tallennettava data
     nc_fname: tiedostonimi
    '''
    try:
        os.remove(nc_fname)
    except OSError as e:
        pass

    # Loop through nowcast timesteps and copy metadata
    # from input POT grib file
    with GribFile(image_grib_file) as grib:
        print("grib:",grib)
        for msg in grib:
            print("msg:",msg)
            for i in range(deterministic_forecast.shape[0]):
                timestep=PD["run_options"]["nowcast_timestep"]
                valid_time = startdate + (i + 1) * dt.timedelta(minutes=timestep)
                print("valid_time:",valid_time)
                msg["bitsPerValue"] = 24
                msg["dataDate"] = int(dt.datetime.strftime(valid_time, "%Y%m%d"))
                msg["dataTime"] = int(dt.datetime.strftime(valid_time, "%H%M"))
                msg["generatingProcessIdentifier"] = 202
                msg["centre"] = 86
                msg["bitmapPresent"] = True
                msg["values"] = deterministic_forecast[i,:,:].flatten()
                with open(str(nc_fname), "ab") as out:
                    msg.write(out)
                    

def calc_datatime(i,nowcast_timestep):
    tdiff_mins = i*int(nowcast_timestep)
    tdiff_hours = int(tdiff_mins/60)
    tdiff_mins = tdiff_mins % (tdiff_hours*60)
    return f"{tdiff_hours:2d}{tdiff_mins:2d}"
                    
def generate_deterministic(observations, motion_field, nowcaster, nowcast_kwargs=None):
    """Generate a deterministic nowcast using semilagrangian extrapolation"""
    # Extrapolation scheme doesn't use the same nowcast_kwargs as steps
    if nowcast_kwargs is None:
        nowcast_kwargs = dict()
        
    forecast = nowcaster(observations, motion_field, PD["run_options"]["leadtimes"],
                         **nowcast_kwargs)
    
    return forecast


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
