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

import numpy as np
from netCDF4 import num2date, Dataset
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt


def arome_model_wind(infile, outdir,
                     vmin=0., vmax=25.4, vmin_pal=0., vmax_pal=50*0.514):
    """
    """
    # Read/Process data
    windfield = Dataset(infile)
    time = windfield.variables['time'][:]
    time_units = windfield.variables['time'].units
    lon = windfield.variables['longitude'][:]
    lat = windfield.variables['latitude'][:]
    #lon0 = lon[0]
    #lat0 = lat[0]
    #dlon = lon[1]-lon[0]
    #dlat = lat[1]-lat[0]
    lon0 = -8.
    lat0 = 53.
    dlon = 0.025
    dlat = -0.025
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    for itime in range(time.size):
        dtime = num2date(time[itime], time_units)
        u10 = windfield.variables['u10'][itime, :, :]
        v10 = windfield.variables['v10'][itime, :, :]
        metadata = {}
        metadata['product_name'] = 'AROME_model_wind'
        metadata['name'] = 'AROME_'+dtime.strftime('%Y%m%dT%HZ')
        metadata['datetime'] = stfmt.format_time(dtime)
        metadata['time_range'] = ['-90m', '+90m']
        metadata['source_URI'] = infile
        metadata['source_provider'] = 'METEO FRANCE'
        metadata['processing_center'] = ''
        metadata['conversion_software'] = 'Syntool'
        metadata['conversion_version'] = '0.0.0'
        metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
        metadata['parameter'] = ['wind speed', 'wind direction']
        metadata['type'] = 'model'
        metadata['model_longitude_resolution'] = abs(dlon)
        metadata['model_latitude_resolution'] = abs(dlat)
        #metadata['model_analysis_datetime'] = stfmt.format_time(rundtime)
        geolocation = {}
        geolocation['projection'] = stfmt.format_gdalprojection()
        geolocation['geotransform'] = [lon[0]-dlon/2., dlon, 0,
                                       lat[0]-dlat/2., 0, dlat]
        band = []
        indndv = np.where(np.ma.getmaskarray(u10) | np.ma.getmaskarray(v10))
        windspeed = np.sqrt(u10**2 + v10**2)
        winddirection = np.mod(np.arctan2(v10, u10)*180./np.pi+360., 360.)
        offset, scale = vmin, (vmax-vmin)/254.
        np.clip(windspeed, vmin, vmax, out=windspeed)
        array = np.round((windspeed - offset) / scale).astype('uint8')
        array[indndv] = 255
        colortable = stfmt.format_colortable('noaa_wind', vmax=vmax, vmax_pal=vmax_pal,
                                             vmin=vmin, vmin_pal=vmin_pal)
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'wind speed', 'unittype':'m/s',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax],
                     'colortable':colortable})
        array = np.round(winddirection/360.*254.).astype('uint8')
        array[indndv] = 255
        band.append({'array':array, 'scale':360./254., 'offset':0.,
                     'description':'wind direction', 'unittype':'deg',
                     'nodatavalue':255, 'parameter_range':[0, 360.]})
        # Write geotiff
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)

