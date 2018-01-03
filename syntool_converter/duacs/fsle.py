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


L4_MAPS = {'FSLE':
           {'productname': 'FSLE',
            'hname': 'fsle_max',
            'parameter': 'Finite Size Lyapunov Exponents',
            'timerange': ['-48h', '+48h']
            },
           }


def fsle_gridded(infile, outdir,
                 vmin=-1., vmax=0., vmin_pal=-1., vmax_pal=0.,
                 write_netcdf=False):
    """
    """
    # Read/Process data
    print('Read/Process data')
    ncfile = NCFile(infile)
    if 'id' in ncfile.read_global_attributes():
        l4id = ncfile.read_global_attribute('id')
    elif (re.match(r'^dt_global_allsat_madt_fsle.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'FSLE'
    else:
        raise Exception('Unknown file.')
    h = ncfile.read_values(L4_MAPS[l4id]['hname'])[0, ::-1, :]
    lon = ncfile.read_values('lon')[:].astype('float64')
    lat = ncfile.read_values('lat')[-1:-3:-1].astype('float64')
    lon[lon > 180.] = lon[lon > 180.]-360.
    indsorted = np.argsort(lon)
    lon = lon[indsorted]
    h = h[:, indsorted]

    for i in range(2):  # avoid rounding errors
        lon[i] = np.round(lon[i]*10000)/10000
        lat[i] = np.round(lat[i]*10000)/10000
    lon0, dlon, lat0, dlat = lon[0], lon[1]-lon[0], lat[0], lat[1]-lat[0]
    dlon = 0.04
    dlat = -0.04
    dtime_units = ncfile.read_field('time').units
    dtime = num2date(ncfile.read_values('time')[0], dtime_units)
    # Construct metadata/geolocation/band(s)
    print('Construct metadata/geolocation/band(s)')
    metadata = {}
    metadata['product_name'] = L4_MAPS[l4id]['productname']
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = L4_MAPS[l4id]['timerange']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'aviso@altimetry.fr'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = [L4_MAPS[l4id]['parameter']]
    # metadata['type'] = 'model'
    # metadata['model_longitude_resolution'] = abs(dlon)
    # metadata['model_latitude_resolution'] = abs(dlat)
    # metadata['model_analysis_datetime'] = stfmt.format_time(rundtime)
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon0-dlon/2., dlon, 0,
                                   lat0-dlat/2., 0, dlat]
    band = []
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(h, vmin, vmax, out=h)
    array = np.round((h - offset) / scale).astype('uint8')
    array[h.mask] = 255
    colortable = stfmt.format_colortable('matplotlib_jet_r',
                                         vmin=vmin, vmax=vmax,
                                         vmin_pal=vmin_pal, vmax_pal=vmax_pal)
    band.append({'array': array, 'scale': scale, 'offset': offset,
                 'description': 'FSLE', 'unittype': 'day',
                 'nodatavalue': 255, 'parameter_range': [vmin, vmax],
                 'colortable': colortable})
    # Write geotiff
    print('Write geotiff')
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
