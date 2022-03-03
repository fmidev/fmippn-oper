"""Utility functions for FMI-PPN"""
import datetime as dt
import numpy as np
import h5py

def pack_value(original_value, scale_factor, add_offset):
    # scale_factor == gain, add_offset == offset
    packed = (original_value - add_offset) / scale_factor
    return packed

def unpack_value(packed_value, scale_factor, add_offset):
    # scale_factor == gain == 1./scaler, add_offset == offset == scale_zero
    original = scale_factor * packed_value + add_offset
    return original

def quantity_is_dbzh(quantity):
    return quantity in {"DBZH", "dbz"}

def quantity_is_rate(quantity):
    return quantity in {"RATE", "rrate"}

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

def prepare_data_for_writing(forecast, options, forecast_undetect=None, forecast_nodata=None):
    """Scale and convert nowcast data to correct datatype.

    The data will be scaled according to equation
        scaled_forecast = (forecast - offset) / gain

    For non-float store datatypes, NaN values are converted to `store_nodata_value`.

    Args:
        forecast (numpy.array): Nowcast data
        options (dict): output options from config

    Returns:
        tuple: (scaled_data, metadata) or (None, dict()) if forecast is None
    """
    # If no forecast, then no need to do anything
    if forecast is None:
        return None, dict()

    # Store data in other datatype format to save space (e.g. float64 -> uint16)
    # If no dtype is given, then default to not converting
    store_dtype = options.get('convert_to_dtype', None)
    store_dtype = np.dtype(store_dtype) if store_dtype is not None else forecast.dtype

    gain = options.get('gain', None)
    if gain is None and options.get("scaler", 0) != 0:
        gain = 1./options["scaler"]

    scale_zero = options.get('offset') if 'offset' in options else options.get("scale_zero")
    if scale_zero in {None, "auto"}:
        scale_zero = np.nanmin(forecast)

    if forecast_nodata is not None:
        store_nodata_value = forecast_nodata
    else:
        cfg_nodata = options.get('set_nodata_value_to', "default")
        store_nodata_value = _get_default_nodata(store_dtype, cfg_nodata)
        
    # Undetect value from input is used in thresholding the data, so let's store that if provided
    # see generate() in ppn.py
    if forecast_undetect is not None:
        undetect = pack_value(forecast_undetect, gain, scale_zero)
    else:
        undetect = options.get('set_undetect_value_to')
        
    # TODO: Take actual forecast_nodata value into account here, if provided
    nodata_mask = ~np.isfinite(forecast)
    fct_scaled = pack_value(forecast, scale_factor=gain, add_offset=scale_zero)
    fct_scaled[nodata_mask] = store_nodata_value
    
    # Tuuli added masking and filling undetect value:
    undetect_mask = (fct_scaled == undetect)
    fct_scaled[undetect_mask] = options.get('set_undetect_value_to')
    undetect = options.get('set_undetect_value_to')
    
    prepared_forecast = fct_scaled.astype(store_dtype)

    metadata = {
        "nodata": store_nodata_value,
        "gain": gain,
        "offset": scale_zero,
        "undetect": undetect,
    }

    return prepared_forecast, metadata

def _get_default_nodata(store_dtype, cfg_nodata):
    if isinstance(cfg_nodata, (int, float)):
        return cfg_nodata

    # If config says "default", then pick one option based on dtype
    default_values = {
        "u": "max_int",  # unsigned int -> maximum value is safer than 0
        "i": "min_int",  # signed int -> minimum value is safer choice than maximum
        "f": "nan"       # float -> nan is available as a special value
    }

    if cfg_nodata == "default":
        cfg_nodata = default_values.get(store_dtype.kind)

    # Actual options
    if cfg_nodata == "max_int":
        return np.iinfo(store_dtype).max

    if cfg_nodata == "min_int":
        return np.iinfo(store_dtype).min

    if cfg_nodata.lower() == "nan":
        return np.nan

    # For unknown or invalid options, raise
    msg = ("Invalid nodata value '{}'. It must be one of following: an integer, a float, 'default',"
           " 'max_int', 'min_int', or 'nan'.".format(cfg_nodata))
    raise TypeError(msg)

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

    # h5py<=2.10 can store variable-length strings in attrs as bytes,
    # while h5py>=3.0 stores them only as strings (fixed-length strings are stored as bytes)
    # ODIM format specifies that they should be bytes, however not every file
    quantity_as_bytes = quantity.encode()
    for attrs in attrs_list:
        if attrs.get('quantity', b'') in {quantity_as_bytes, quantity}:
            und = attrs.get('undetect', None)
            if und is None:
                raise RuntimeError(f"'undetect' attribute is missing from {quantity} data attributes!")
            gain = attrs.get('gain', 1.0)
            offset = attrs.get('offset', 0.0)
            return unpack_value(und, gain, offset)

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

    # Deduce ODIM format quantity from units
    unit = metadata.get("unit", "Unknown").lower()
    if unit == "dbz":
        quantity="DBZH"
    elif unit in {"rrate", "mm/h"}:
        quantity="RATE"
    else:
        quantity = "Unknown"

    #Create data/what group and store metadata
    data_what_grp=data_grp.create_group("what")
    data_what_grp.attrs["quantity"] = quantity
    data_what_grp.attrs["gain"] = scale_meta.get("gain")
    data_what_grp.attrs["offset"] = scale_meta.get("offset")
    data_what_grp.attrs["nodata"] = scale_meta.get("nodata")
    data_what_grp.attrs["undetect"] = scale_meta.get("undetect")
