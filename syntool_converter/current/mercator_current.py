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
import os
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime


# infile = '/mnt/data/mercator/global-analysis-forecast-phys-001-002_agulhas_daily.nc'
# outdir = '/mnt/data/syntool_inputs/ovl'
# date0 = datetime(2015,3,1,12)
# for i in range(31): mercator_current(infile, outdir, date=date0 + i * timedelta(hours=24))

def mercator_current(infile, outdir, date=None,
                     vmin=0., vmax=5.08, vmin_pal=0., vmax_pal=2.):
    """ """
    if date is None:
        raise Exception('mercator_current conversion needs a date !')
    # Read/Process data
    print 'Read/Process data'
    ncfile = Dataset(infile)
    time = ncfile.variables['time_counter']
    time_index = np.where(num2date(time[:], time.units) == date)[0]
    if time_index.size != 1:
        raise Exception('Date not found in mercator file !')
    time_index = time_index[0]
    dtime = num2date(time[time_index], time.units)
    lat = ncfile.variables['latitude'][::-1]
    lon = ncfile.variables['longitude'][:]
    ucur = ncfile.variables['u'][time_index, 0, ::-1, :]
    vcur = ncfile.variables['v'][time_index, 0, ::-1, :]
    if isinstance(ucur, np.ma.MaskedArray) == False:
        ucur = np.ma.masked_invalid(ucur)
    if isinstance(vcur, np.ma.MaskedArray) == False:
        vcur = np.ma.masked_invalid(vcur)
    if 'daily' in os.path.basename(infile) and \
       'agulhas' in os.path.basename(infile):
        product_name = 'MERCATOR_current_daily_agulhas'
        time_range = ['-12h', '+12h']
    else:
        raise Exception('Mercator product not taken into account for now.')
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = product_name
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0] + '_' +\
                       dtime.strftime('%Y%m%dT%H%M%S')
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'MERCATOR'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = ['current velocity', 'current direction']
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    geolocation['geotransform'] = [lon[0]-dlon/2., dlon, 0,
                                   lat[0]-dlat/2., 0, dlat]
    band = []
    mask = np.ma.getmaskarray(ucur) | np.ma.getmaskarray(vcur)
    curvel = np.sqrt(ucur.data**2 + vcur.data**2)
    curdir = np.mod(np.arctan2(vcur.data, ucur.data)*180./np.pi+360., 360.)
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
