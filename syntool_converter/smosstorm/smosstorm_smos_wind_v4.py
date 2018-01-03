#!/usr/bin/env python
# coding=utf-8

"""
Copyright (C) 2014-2018 OceanDataLab

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from netCDF4 import Dataset, num2date
import numpy as np
from datetime import datetime, timedelta
import syntool_converter.utils.syntoolformat as stfmt
import os


# def check_flags():
#     """
#     """
#     import matplotlib.pyplot as plt
#     import glob
#     paths = glob.glob('/mnt/data/smosstorm/data/smosstorm/l2/smos/*/*/*.nc')
#     for path in paths:
#         dataset = Dataset(path)
#         # print len(dataset.dimensions['lon']), os.path.basename(path)
#         # dataset.close()
#         # continue
#         wind = dataset.variables['wind_speed'][0, :, :]
#         flags = dataset.variables['flags'][0, :, :]
#         flags_bytes = [0, 3, 4, 5, 6] # advised in product
#         flags_mask = np.any([np.bitwise_and(flags, 2 ** b) != 0 for b in flags_bytes], axis=0)
#         mask = np.ma.getmaskarray(wind) | flags_mask
#         wind = np.ma.MaskedArray(np.ma.getdata(wind), mask=mask)
#         vmin, vmax = wind.compressed().min(), wind.compressed().max()
#         plt.figure(figsize=(16., 12.))
#         plt.subplot(1, 2, 1)
#         plt.imshow(dataset.variables['wind_speed'][0, :, :], interpolation='nearest',
#                    vmin=vmin, vmax=vmax)
#         plt.colorbar()
#         plt.subplot(1, 2, 2)
#         plt.imshow(wind, interpolation='nearest', vmin=vmin, vmax=vmax)
#         plt.colorbar()
#         dataset.close()
#         plt.show()
#         import pdb ; pdb.set_trace()


def smosstorm_smos_wind_v4(infile, outdir,
                           vmin=0., vmax=50.8, vmin_pal=0., vmax_pal=50.):
    """
    """
    # Read/Process data
    dataset = Dataset(infile)
    if len(dataset.dimensions['time']) == 0:
        raise Exception('time dimension of null length.')
    # Get masked wind using appropriate flags
    wind = dataset.variables['wind_speed'][0, :, :]
    flags = dataset.variables['flags'][0, :, :]
    flags_bytes = [0, 3, 4, 5, 6] # advised in product
    flags_mask = np.any([np.bitwise_and(flags, 2 ** b) != 0 for b in flags_bytes], axis=0)
    mask = np.ma.getmaskarray(wind) | flags_mask
    if mask.all():
        raise Exception('Data is fully masked.')
    wind = np.ma.MaskedArray(np.ma.getdata(wind), mask=mask)
    # Get lon/lat
    lat = dataset.variables['lat'][:]
    dlat = lat[1] - lat[0]
    if not np.all(np.isclose(lat[1:] - lat[:-1], dlat)):
        raise Exception('Unexpected not unique dlat.')
    lon = dataset.variables['lon'][:]
    dlon = lon[1] - lon[0]
    if not np.all(np.isclose(lon[1:] - lon[:-1], dlon)):
        raise Exception('Unexpected not unique dlon.')
    # Get time
    time_units = dataset.variables['time'].units
    if 'dtime' in dataset.variables:
        # standard case
        if 'days' not in time_units:
            raise Exception('time units expected in days.')
        if 'days' not in dataset.variables['dtime'].units:
            raise Exception('dtime units expected in days.')
        time = dataset.variables['time'][0] + dataset.variables['dtime'][0, :, :]
    else:
        # rare case
        time = dataset.variables['time'][0, :, :]
    _start_time = datetime.strptime(dataset.time_coverage_start, '%Y%m%dT%H%M%SZ')
    _stop_time = datetime.strptime(dataset.time_coverage_end, '%Y%m%dT%H%M%SZ')
    dataset.close()
    # Keep only valid part of grid
    valid = np.where(np.ma.getmaskarray(wind) == False)
    lat_slice = slice(valid[0].min(), valid[0].max() + 1)
    lon_slice = slice(valid[1].min(), valid[1].max() + 1)
    wind = wind[lat_slice, lon_slice]
    lat = lat[lat_slice]
    lon = lon[lon_slice]
    time = time[lat_slice, lon_slice]
    # Set start_time/stop_time
    if np.ma.is_masked(time) and time.mask.all():
        raise Exception('time is fully masked in valid slice.')
    start_time = num2date(time.min(), time_units)
    stop_time = num2date(time.max(), time_units)
    del time
    if (start_time + timedelta(seconds=1)) < _start_time or \
       (stop_time - timedelta(seconds=1)) > _stop_time:
        raise Exception('time outside time coverage in global attributes.')
    # Flip grid
    if dlat > 0:
        wind = wind[::-1, :]
        lat = lat[::-1]
        dlat = -dlat
    if dlon < 0:
        wind = wind[:, ::-1]
        lon = lon[::-1]
        dlon = -dlon
    # Rearrange grid
    # (shift the left part of the grid to minimize valid extent)
    lon_min, lon_max = lon[0], lon[-1]
    if lon_min < -180. or lon_max > 180.:
        raise Exception('Unexpected lon outside [-180, 180].')
    valid_lon = lon[np.any(~np.ma.getmaskarray(wind), axis=0)]
    if valid_lon.size == 1:
        shifted_lon_min = lon_min - 1
        shifted_lon_max = lon_max + 1
    else:
        valid_dlon = valid_lon[1:] - valid_lon[:-1]
        shifted_lon_min = valid_lon[valid_dlon.argmax() + 1]
        shifted_lon_max = valid_lon[valid_dlon.argmax()] + 360.
    if (lon_max - lon_min) > (shifted_lon_max - shifted_lon_min):
        _wind = wind.copy()
        shifted_nlon = np.round((shifted_lon_max - shifted_lon_min) / dlon).astype('int') + 1
        wind = np.ma.masked_all((_wind.shape[0], shifted_nlon), dtype=_wind.dtype)
        indlon = np.round(np.mod(lon - shifted_lon_min, 360.) / dlon).astype('int')
        ind = np.where(indlon < shifted_nlon)[0]
        wind[:, indlon[ind]] = _wind[:, ind]
        lon = np.linspace(shifted_lon_min, shifted_lon_max, num=shifted_nlon, endpoint=True)

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SMOSSTORM_SMOS_wind_V4'
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'ESA'
    metadata['processing_center'] = 'IFREMER'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'wind speed'
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon[0] - dlon / 2., dlon, 0,
                                   lat[0] - dlat / 2., 0, dlat]
    band = []
    indndv = np.where(np.ma.getmaskarray(wind))
    offset, scale = vmin, (vmax - vmin) / 254.
    np.clip(wind.data, vmin, vmax, out=wind.data)
    array = np.round((wind.data - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('matplotlib_jet',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'wind speed', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
