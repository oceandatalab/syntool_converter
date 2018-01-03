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

from cerbere.mapper.ghrsstncfile import GHRSSTNCFile
import numpy as np
from netCDF4 import num2date
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os


def ghrsst_seviri(infile, outdir,
                  vmin=271.05, vmax=309.15, vmin_pal=273., vmax_pal=305.,
                  write_netcdf=False):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    seviri = GHRSSTNCFile(infile)
    sst = seviri.read_values('sea_surface_temperature')[::-1, :]
    mask = seviri.read_values('quality_level')[::-1, :]
    #sea_ice_fraction = seviri.read_values('sea_ice_fraction')[::-1, :]
    # lon = seviri.read_values('lon')
    # dlon = lon[1] - lon[0]
    # lon0 = lon[0] - dlon / 2
    # lat = seviri.read_values('lat')[::-1]
    # dlat = lat[1] - lat[0]
    # lat0 = lat[0] - dlat / 2
    lon0 = seviri.read_global_attribute('westernmost_longitude')
    dlon = float(seviri.read_global_attribute('geospatial_lon_resolution'))
    lat0 = seviri.read_global_attribute('northernmost_latitude')
    dlat = -float(seviri.read_global_attribute('geospatial_lat_resolution'))
    dtime = seviri.read_values('time')[0]
    dtime_units = seviri.read_field('time').units
    dtime = num2date(dtime, dtime_units)
    start_time = datetime.strptime(seviri.read_global_attribute('start_time'),
                                   '%Y%m%dT%H%M%SZ')
    stop_time = datetime.strptime(seviri.read_global_attribute('stop_time'),
                                  '%Y%m%dT%H%M%SZ')
    (dtime, time_range) = stfmt.format_time_and_range(start_time,
                                                      stop_time,
                                                      units='ms')
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = 'SEVIRI_SST'
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'Ifremer'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea surface temperature'
    metadata['type'] = 'remote sensing'
    metadata['longitude_resolution'] = abs(dlon)
    metadata['latitude_resolution'] = abs(dlat)
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon0, dlon, 0,
                                   lat0, 0, dlat]
    band = []
    #indndv = np.where((sst.mask == True) | (sea_ice_fraction > 0))
    indndv = np.where((sst.mask == True) | (mask <= 3))
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(sst, vmin, vmax, out=sst)
    array = np.round((sst - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('cerbere_medspiration',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea surface temperature', 'unittype':'K',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
    if write_netcdf == False:
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    elif write_netcdf == True:
        print 'Write netcdf'
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        band[0]['name'] = 'sea_surface_temperature'
        band[0]['long_name'] = 'sea surface subskin temperature'
        band[0]['standard_name'] = 'sea_surface_subskin_temperature'
        metadata['spatial_resolution'] = min([abs(dlat), abs(dlon)]) * 111000.
        dgcps = np.round(1. / np.abs(np.array([dlat, dlon]))).astype('int')
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'grid_lonlat',
                           dgcps=dgcps)

