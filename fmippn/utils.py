"""Utility functions for FMI-PPN"""
import datetime as dt
import numpy as np
import h5py


def utcnow_floored(increment=5):
    """Return UTC time with minutes replaced by latest multiple of `increment`."""
    now = dt.datetime.utcnow()
    floored_minutes = now.minute - (now.minute % increment)
    now = now.replace(minute=floored_minutes)
    return now


def store_timeseries(grp, data, startdate, timestep, metadata=None):
    """Store timeseries for one nowcast ensemble member.

    Input:
        grp -- group object in HDF5 file
        data -- timeseries data in 3-dimensional numpy array (first dimension=time)
        startdate -- nowcast analysis time (datetime object)
        timestep -- time difference between nowcast fields (int)

    Optional input:
        metadata -- a dictionary containing additional metadata. Will be added
                    to dataset as attributes
    """
    if metadata is None:
        metadata = dict()
    for index in range(data.shape[0]):
        ts_point = data[index, :, :]
        tmp = grp.create_dataset("leadtime-{:0>2}".format(index), data=ts_point)
        valid_time = startdate + (index + 1) * dt.timedelta(minutes=timestep)
        tmp.attrs["Valid for"] = int(dt.datetime.strftime(valid_time, "%Y%m%d%H%M%S"))
        for key, value in metadata.items():
            tmp.attrs[key] = value


def prepare_fct_for_saving(fct, scaler, scale_zero, store_dtype, store_nodata_value):
    """Scale and convert nowcast data to correct datatype.

    The data will be scaled according to equation
        fct_scaled = scaler * (fct - scale_zero)

    NaN values are converted to `store_nodata_value`.

    Input:
        fct -- nowcast data (numpy.array)
        scaler -- scaling term for scale equation
        scale_zero -- value for data that will be 0 in scaled values
        store_dtype -- data type for scaled values (numpy.dtype)
        store_nodata_value -- Value that will be used to mark invalid elements
                              in scaled values

    Output:
        fct_scaled -- scaled and converted data
    """
    nodata_mask = ~np.isfinite(fct)
    fct_scaled = scaler * (fct - scale_zero)
    if store_nodata_value != -1 and np.any(fct_scaled >= store_nodata_value):
        raise ValueError("Cannot store forecast to a file: One or more values would be "
                         "larger than maximum allowed value (%i) causing overflow. "
                         % (store_nodata_value-1))
    fct_scaled[nodata_mask] = store_nodata_value
    fct_scaled = fct_scaled.astype(store_dtype)
    return fct_scaled


def copy_odim_attributes(infile,outf):
    """Copy attribute groups /what, /where and /how from
    input ODIM HDF5 file to output ODIM HDF5 file as they are.
    
    Keyword arguments:
    infile -- ODIM HDF5 input composite filename
    outf -- FMI-PPN output HDF5 file object
    """

    inf=h5py.File(infile, 'r')

    #Copy attribute groups /what, /where and /how
    what=outf.create_group("what")
    for key, val in inf["what"].attrs.items():
         what.attrs[key] = val

    where=outf.create_group("where")
    for key, val in inf["where"].attrs.items():
        where.attrs[key] = val

    how=outf.create_group("how")
    for key, val in inf["how"].attrs.items():
        how.attrs[key] = val
        
    inf.close()



    
