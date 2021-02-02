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

def get_odim_attrs_from_input(infile):
    """Read attribute groups /what, /where and /how from input ODIM HDF5 file.

    Input:
        infile -- ODIM HDF5 input composite filename

    Output:
        A dictionary with dictionaries 'what', 'where' and 'how' containing the attributes.
    """
    with h5py.File(infile, 'r') as f:
        what = dict(f["what"].attrs)
        where = dict(f["where"].attrs)
        how = dict(f["how"].attrs)

    return {
        "what": what,
        "where": where,
        "how": how,
    }

def get_odim_data_undetect(fname, quantity='DBZH'):
    """

    Input:
        fname --
        quantity -- (default: DBZH)

    Output:
        Value for 'undetect' pixels in actual units (using formula 'gain * undetect + offset')

    Raises:
        RuntimeError if value for 'undetect' cannot be found or is missing from attributes
    """
    # This assumes that there is only one 'dataset/data' group containing the wanted quantity
    # TODO: Generalize if there can be more groups with the same quantity
    attrs_list = []
    def _append_attrs(name, obj):
        """Append attribute dictionary to attrs_list (this abuses a SIDE EFFECT!)"""
        if 'quantity' in obj.attrs:
            attrs_list.append(dict(obj.attrs))
        return None

    # h5py.visititems goes through the file recursively
    # it needs a callable with two arguments (name, obj)
    # See h5py documentation for more
    with h5py.File(fname, 'r') as f:
        f.visititems(_append_attrs)

    for attrs in attrs_list:
        if attrs.get('quantity', b'').decode() == quantity:
            und = attrs.get('undetect', None)
            if und is None:
                raise RuntimeError(f"'undetect' attribute is missing from {quantity} data attributes!")
            gain = attrs.get('gain', 1.0)
            offset = attrs.get('offset', 0.0)
            return gain * und + offset

    raise RuntimeError(f"Could not find 'undetect' value for {quantity} in file {fname}")


def copy_odim_attributes(odim_metadata,outf):
    """Copy attribute groups /what, /where and /how from
    input ODIM HDF5 file to output ODIM HDF5 file as they are.

    Keyword arguments:
    odim_metadata -- dictionary containing subdictionaries what, where 
                     and how read from ODIM HDF5 input composite
    outf -- FMI-PPN output HDF5 file object
    """
    
    #Copy attribute groups /what, /where and /how
    what=outf.create_group("what")
    for key, val in odim_metadata["what"].items():
        what.attrs[key] = val

    where=outf.create_group("where")
    for key, val in odim_metadata["where"].items():
        where.attrs[key] = val

    how=outf.create_group("how")
    for key, val in odim_metadata["how"].items():
        how.attrs[key] = val
   
    

def store_odim_dset_attrs(dset_grp, dset_index, startdate, timestep):
    """Store ODIM attributes to datasets. Each dataset
    represents a different timestep.

    Keyword arguments:
    dset_grp -- dataset HDF5 group object
    dset_index -- dataset number (dataset1, dataset2 etc), indicates number of timestep
    startdate -- nowcast analysis time (datetime object)
    timestep -- time difference between nowcast fields (int)
    """

    #Calculate valid time for each step
    valid_time = startdate + (dset_index + 1) * dt.timedelta(minutes=timestep)

    #Add attributes to each dataset
    dset_how_grp=dset_grp.create_group("how")
    dset_how_grp.attrs["simulated"]="True"

    dset_what_grp=dset_grp.create_group("what")
    dset_what_grp.attrs["startdate"] = int(dt.datetime.strftime(valid_time, "%Y%m%d"))
    dset_what_grp.attrs["enddate"] = int(dt.datetime.strftime(valid_time, "%Y%m%d"))
    dset_what_grp.attrs["starttime"] = int(dt.datetime.strftime(valid_time, "%H%M%S"))
    dset_what_grp.attrs["endtime"] = int(dt.datetime.strftime(valid_time, "%H%M%S"))


def store_odim_data_what_attrs(data_grp,metadata,scale_meta):
    """Store ODIM attributes to data/what group. Each data group (data1, data2 ...)
    represents a different ensemble member.

    Keyword arguments:
    data_grp -- data HDF5 group object (data1, data2 ...)
    metadata -- array containing FMIPPN output metadata
    scale_meta -- scale values metadata
    """

    #Change quantity to ODIM format
    quantity=metadata.get("unit","Unknown")
    if quantity == "dBZ":
        quantity="DBZH"
    elif quantity == "rrate":
        quantity="RATE"

    #Create data/what group and store metadata
    data_what_grp=data_grp.create_group("what")
    data_what_grp.attrs["quantity"] = quantity
    data_what_grp.attrs["gain"] = scale_meta.get("gain")
    data_what_grp.attrs["offset"] = scale_meta.get("offset")
    data_what_grp.attrs["nodata"] = scale_meta.get("nodata")
    data_what_grp.attrs["undetect"] = scale_meta.get("undetect")
