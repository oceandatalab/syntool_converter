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


L4_MAPS = {'Surface_height':
           {'productname': 'Surface_height_current',
            'uname': 'u',
            'vname': 'v',
            'parameter': 'Current from sea surface height',
            'timerange': ['0s', '+86399s']
            },
           'Sea_Level_Anomaly':
           {'productname': 'Sea_Level_Anomaly_current',
            'uname': 'u',
            'vname': 'v',
            'parameter': 'Current from Sea level anomaly',
            'timerange': ['0s', '+86399s']
            },
           'Mean_Dynamic_Topo':
           {'productname': 'Mean_Dynamic_Topography_current',
            'uname': 'u',
            'vname': 'v',
            'parameter': 'Current from Mean Dynamic Topography',
            'timerange': ['-31d', '+31d']
            },
           'Tide':
           {'productname': 'barotropic_tidal_height_current',
            'uname': 'u',
            'vname': 'v',
            'parameter': 'Current tides',
            'timerange': ['-30m', '+30m']
            },
           }


def current(infile, outdir,
            vmin=0., vmax=1.50, vmin_pal=0., vmax_pal=1.5,
            write_netcdf=False):
    """
    """
    # Read/Process data
    print('Read/Process data')
    ncfile = NCFile(infile)
    if 'id' in ncfile.read_global_attributes():
        l4id = ncfile.read_global_attribute('id')
    elif (re.match(r'^dt_global_allsat_madt_uv.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Surface_height'
        # vmin = 0.; vmax = 2.; vmin_pal = 0.; vmax_pal = 2.
    elif (re.match(r'^dt_global_allsat_msla_uv.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Sea_Level_Anomaly'
        # vmin = 0; vmax = 1; vmin_pal = 0.; vmax_pal = 1
    elif re.match(r'^mdt.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'Mean_Dynamic_Topo'
        # vmin = 0; vmax = 1.5; vmin_pal = 0; vmax_pal = 1.5
    elif re.match(r'^Tide_.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'Tide'
        # vmin = 0; vmax = 1.5; vmin_pal = 0; vmax_pal = 1.5
    else:
        raise Exception('Unknown file.')
    # /TMP
    ucur = ncfile.read_values(L4_MAPS[l4id]['uname'])[0, ::-1, :]
    vcur = ncfile.read_values(L4_MAPS[l4id]['vname'])[0, ::-1, :]
    lon = ncfile.read_values('lon')[:].astype('float64')
    lat = ncfile.read_values('lat')[-1:-3:-1].astype('float64')
    lon[lon > 180.] = lon[lon > 180.]-360.
    indsorted = np.argsort(lon)
    lon = lon[indsorted]
    ucur = ucur[:, indsorted]
    vcur = vcur[:, indsorted]
    lon0, dlon, lat0, dlat = lon[0], lon[1]-lon[0], lat[0], lat[1]-lat[0]
    if l4id in ['Mean_Dynamic_Topo', ]:
        dtime = datetime(2014, 12, 1)
    else:
        dtime_units = ncfile.read_field('time').units
        dtime = num2date(ncfile.read_values('time')[0], dtime_units)
    # rundtime = ncfile.read_global_attribute('date_modified')
    # rundtime = datetime.strptime(rundtime, '%Y%m%dT%H%M%SZ')
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
    if l4id in ['Tide', ]:
        ucur = ucur[:-1, ::]
        vcur = vcur[:-1, ::]
        geolocation['geotransform'] = [lon0, dlon, 0, lat0, 0, dlat]
    else:
        geolocation['geotransform'] = [lon0-dlon/2., dlon, 0,
                                       lat0-dlat/2., 0, dlat]
    band = []
    offset, scale = vmin, (vmax-vmin)/254.
    mask = ucur.mask | vcur.mask
    curvel = np.sqrt(ucur.data**2 + vcur.data**2)
    curdir = np.mod(np.arctan2(vcur.data, ucur.data)*180./np.pi+360., 360.)
    np.clip(curvel, vmin, vmax, out=curvel)
    array = np.round((curvel - offset) / scale).astype('uint8')
    array[mask] = 255
    colortable = stfmt.format_colortable('matplotlib_jet', vmin=vmin,
                                         vmax=vmax, vmin_pal=vmin_pal,
                                         vmax_pal=vmax_pal)
    band.append({'array': array, 'scale': scale, 'offset': offset,
                 'description': L4_MAPS[l4id]['productname'],
                 'unittype': 'm/s', 'nodatavalue': 255,
                 'parameter_range': [vmin, vmax], 'colortable': colortable})
    array = np.round(curdir/360.*254.).astype('uint8')
    array[mask] = 255
    band.append({'array': array, 'scale': 360./254., 'offset': 0.,
                 'description': 'current direction', 'unittype': 'deg',
                 'nodatavalue': 255, 'parameter_range': [0, 360.]})

    # Write geotiff
    print('Write geotiff')
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
