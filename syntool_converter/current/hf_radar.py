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

import netCDF4
import numpy as np
import numpy.ma
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt


# US West Coast
# 'http://hfrnet.ucsd.edu/thredds/dodsC/HFRNet/USWC/6km/hourly/RTV'
# US East and Gulf Coast
# 'http://hfrnet.ucsd.edu/thredds/dodsC/HFRNet/USEGC/6km/hourly/RTV'
# US Hawai'i State
# 'http://hfrnet.ucsd.edu/thredds/dodsC/HFRNet/USHI/6km/hourly/RTV'
# Alaska - North Slope
# 'http://hfrnet.ucsd.edu/thredds/dodsC/HFRNet/AKNS/6km/hourly/RTV'
# Puerto Rico and the US Virgin Islands
# 'http://hfrnet.ucsd.edu/thredds/dodsC/HFRNet/PRVI/6km/hourly/RTV'


# def search_in_time(time, time_search):
#     """
#     """
#     indb = 0
#     inde = len(time) - 1
#     found = False
#     while indb <= inde and found == False:
#         indm = (indb + inde) / 2
#         timem = time[indm]
#         if timem == time_search:
#             found = True
#         elif timem < time_search:
#             indb = indm + 1
#         elif timem > time_search:
#             inde = indm - 1
#     if found == False:
#         indm = None
#     return indm


def hf_radar(url, outdir, date=None,
             vmin=0., vmax=5.08, vmin_pal=0., vmax_pal=2.):
    """
    """
    if date is None:
        raise Exception('A date has to be specified for HF radar !')
    # Read/Process data
    print 'Read/Process data'
    dataset = netCDF4.Dataset(url+'?lat,lon,time,u,v', 'r')
    time = dataset.variables['time']
    time_search = np.round(netCDF4.date2num(date, time.units))
    # time_index = search_in_time(time, time_search)
    # if time_index is None:
    #     raise Exception('Date not found !')
    time_index = np.where((time[:] == time_search))[0]
    if time_index.size == 0:
        raise Exception('Date not found !')
    time_index = time_index[0]
    dtime = netCDF4.num2date(time[time_index], time.units)
    lon = dataset.variables['lon'][:]
    dlon = (lon[-1] - lon[0]) / (lon.size - 1)
    lat = dataset.variables['lat'][::-1]
    dlat = (lat[-1] - lat[0]) / (lat.size - 1)
    uvel = dataset.variables['u'][time_index, ::-1, :]
    if not isinstance(uvel, numpy.ma.MaskedArray):
        uvel = numpy.ma.masked_invalid(uvel)
    vvel = dataset.variables['v'][time_index, ::-1, :]
    if not isinstance(vvel, numpy.ma.MaskedArray):
        vvel = numpy.ma.masked_invalid(vvel)
    name = '_'.join(url.split('/')[-4:-1]) + dtime.strftime('_%Y%m%dT%H%M%S')

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = 'HF_radar'
    metadata['name'] = name
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = ['-30m', '+30m']
    metadata['source_URI'] = url
    metadata['source_provider'] = 'Scripps Institution of Oceanography'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = ['current velocity', 'current direction']
    # metadata['type'] = 'model'
    # metadata['model_longitude_resolution'] = abs(dlon)
    # metadata['model_latitude_resolution'] = abs(dlat)
    # metadata['model_analysis_datetime'] = stfmt.format_time(rundtime)
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon[0]-dlon/2., dlon, 0,
                                   lat[0]-dlat/2., 0, dlat]
    band = []
    mask = uvel.mask | vvel.mask
    curvel = np.sqrt(uvel.data ** 2  + vvel.data ** 2)
    curdir = np.mod(np.arctan2(vvel.data, uvel.data)*180./np.pi+360., 360.)
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(curvel, vmin, vmax, out=curvel)
    array = np.round((curvel - offset) / scale).astype('uint8')
    array[mask] = 255
    colortable = stfmt.format_colortable('matplotlib_jet', vmin=vmin, vmax=vmax,
                                         vmin_pal=vmin_pal, vmax_pal=vmax_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'current velocity', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    array = np.round(curdir/360.*254.).astype('uint8')
    array[mask] = 255
    band.append({'array':array, 'scale':360./254., 'offset':0.,
                 'description':'current direction', 'unittype':'deg',
                 'nodatavalue':255, 'parameter_range':[0, 360.]})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)

