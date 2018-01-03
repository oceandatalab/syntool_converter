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

from cerbere.mapper.ncfile import NCFile
import numpy as np
from netCDF4 import num2date
import os
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
import re


L4_MAPS = {
    # New format v01
    'e-Odyn':
    {
        'productname': 'e-Odyn_current_from_AIS',
        'uname': 'u',
        'vname': 'v',
        'timerange': ['-12h', '+12h']
    },
    'EODYN-VELOCITY-':
    {
        'productname': 'e-Odyn_current_from_AIS',
        'uname': 'u',
        'vname': 'v',
        'timerange': ['-12h', '+12h']
    },
}


def eodyn_current(infile, outdir,
                   vmin=0., vmax=5.08, vmin_pal=0., vmax_pal=2.,
                   write_netcdf=False):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    ncfile = NCFile(infile)
    if 'id' in ncfile.read_global_attributes():
        l4id = ncfile.read_global_attribute('id')
#        l4id = 'e-Odyn' #ncfile.read_global_attribute('id')
    elif re.match(r'^e-Odyn_.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'e-Odyn'
    else:
        raise Exception('Unknown GlobCurrent L4 file.')
    # /TMP
    ucur = ncfile.read_values(L4_MAPS[l4id]['uname'])[::, ::-1, 0]
    ucur = np.transpose(ucur)
    vcur = ncfile.read_values(L4_MAPS[l4id]['vname'])[::, ::-1, 0]
    vcur = np.transpose(vcur)
    masku = [ucur == -9999]
    maskv = [vcur == -9999]
    if l4id not in ['CourantGeostr']:
        lon = ncfile.read_values('lon')[0:2].astype('float64')
        lat = ncfile.read_values('lat')[-1:-3:-1].astype('float64')
        for i in range(2): # avoid rounding errors
            lon[i] = np.round(lon[i]*10000)/10000
            lat[i] = np.round(lat[i]*10000)/10000
    else:
        lon = ncfile.read_values('lon')[:]
        shift = -np.where(lon < 0)[0][0]
        ucur = np.roll(ucur, shift, axis=1)
        vcur = np.roll(vcur, shift, axis=1)
        lon = lon[shift:shift+2]
        lat = ncfile.read_values('lat')[-1:-3:-1]
    lon0, dlon, lat0, dlat = lon[0], lon[1]-lon[0], lat[0], lat[1]-lat[0]
    #dtime_units = ncfile.read_field('time').units
    #dtime = num2date(ncfile.read_values('time')[0], dtime_units)
    timefmt = '%Y-%m-%dT%H:%M:%S.%fZ'
    start_time = datetime.strptime(ncfile.read_global_attribute('time_coverage_start'),
                                   timefmt)
    stop_time = datetime.strptime(ncfile.read_global_attribute('time_coverage_end'),
                                   timefmt)

    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                              units='ms')
    # rundtime = ncfile.read_global_attribute('date_modified')
    # rundtime = datetime.strptime(rundtime, '%Y%m%dT%H%M%SZ')
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = L4_MAPS[l4id]['productname']
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    #metadata['time_range'] = L4_MAPS[l4id]['timerange']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'e-Odyn'
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
    geolocation['geotransform'] = [lon0-dlon/2., dlon, 0,
                                   lat0-dlat/2., 0, dlat]
    band = []
    mask = ucur.mask | vcur.mask
    print(mask)
    curvel = np.sqrt(ucur.data**2 + vcur.data**2)
    curdir = np.mod(np.arctan2(vcur.data, ucur.data)*180./np.pi+360., 360.)
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(curvel, vmin, vmax, out=curvel)
    array = np.round((curvel - offset) / scale).astype('uint8')
    array[mask] = 255
    array[masku] = 255
    array[maskv] = 255
    print(array)
    colortable = stfmt.format_colortable('matplotlib_jet', vmin=vmin, vmax=vmax,
                                         vmin_pal=vmin_pal, vmax_pal=vmax_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'current velocity', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    array = np.round(curdir/360.*254.).astype('uint8')
    array[mask] = 255
    array[masku] = 255
    array[maskv] = 255
    band.append({'array':array, 'scale':360./254., 'offset':0.,
                 'description':'current direction', 'unittype':'deg',
                 'nodatavalue':255, 'parameter_range':[0, 360.]})
    # Write geotiff
    if write_netcdf == False:
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    elif write_netcdf == True:
        print 'Write netcdf'
        # u/v -> bands
        band = []
        mask = ucur.mask | vcur.mask
        vmin = -vmax
        offset, scale = vmin, (vmax - vmin) / 254.
        u = np.clip(ucur.data, vmin, vmax)
        array = np.round((u - offset) / scale).astype('uint8')
        array[mask] = 255
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'current u', 'unittype':'m/s',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax]})
        v = np.clip(vcur.data, vmin, vmax)
        array = np.round((v - offset) / scale).astype('uint8')
        array[mask] = 255
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'current v', 'unittype':'m/s',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax]})
        # Write
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        stfmt.write_netcdf(ncfile, metadata, geolocation, band,
                           dgcpy=1., dgcpx=1.)
