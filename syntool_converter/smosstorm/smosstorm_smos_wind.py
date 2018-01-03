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

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os


def smosstorm_smos_wind(infile, outdir,
                        vmin=0., vmax=50.8, vmin_pal=0., vmax_pal=50.):
    """
    """
    # Read/Process data
    dataset = Dataset(infile)
    wind = dataset.variables['wind_speed'][0, :, :]
    valid = np.where(ma.getmaskarray(wind) == False)
    lat_slice = slice(valid[0].min(), valid[0].max() + 1)
    lon_slice = slice(valid[1].min(), valid[1].max() + 1)
    nlon = dataset.variables['lon'].shape[0]
    if lon_slice.stop - lon_slice.start < nlon / 2:
        wind = wind[lat_slice, lon_slice]
        lat = dataset.variables['lat'][lat_slice]
        lon = dataset.variables['lon'][lon_slice]
    else:
        lon_first = dataset.variables['lon'][0]
        lon_last = dataset.variables['lon'][-1]
        validl = np.where(ma.getmaskarray(wind[lat_slice, 0:nlon / 2]) == False)
        lon_slicel = slice(0, validl[1].max() + 1)
        validr = np.where(ma.getmaskarray(wind[lat_slice, nlon / 2:]) == False)
        if (lon_last - lon_first) % 360 == 0:
            lon_slicer = slice(validr[1].min() + nlon / 2, nlon - 1)
        else:
            lon_slicer = slice(validr[1].min() + nlon / 2, nlon)
        wind = np.ma.hstack((wind[lat_slice, lon_slicer],
                             wind[lat_slice, lon_slicel]))
        lat = dataset.variables['lat'][lat_slice]
        lon = np.hstack((dataset.variables['lon'][lon_slicer],
                         dataset.variables['lon'][lon_slicel] + 360))
    start_time = datetime.strptime(dataset.time_coverage_start, '%Y%m%dT%H%M%S')
    stop_time = datetime.strptime(dataset.time_coverage_stop, '%Y%m%dT%H%M%S')

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SMOSSTORM_SMOS_wind'
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
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    geolocation['geotransform'] = [lon[0] - dlon / 2., dlon, 0,
                                   lat[0] - dlat / 2., 0, dlat]
    band = []
    indndv = np.where(ma.getmaskarray(wind))
    offset, scale = vmin, (vmax-vmin)/254.
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
