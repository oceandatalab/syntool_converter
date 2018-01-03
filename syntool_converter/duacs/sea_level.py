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


L4_MAPS = {'Mean_Sea_Level':
           {'productname': 'Mean_Sea_Level',
            'hname': 'sea_level_trends',
            'parameter': 'mean sea level',
            'timerange': ['-31d', '+31d']
            },
           'Surface_height':
           {'productname': 'Surface_height',
            'hname': 'adt',
            'parameter': 'surface height',
            'timerange': ['0s', '+86399s']
            },
           'Sea_Level_Anomaly':
           {'productname': 'Sea_Level_Anomaly',
            'hname': 'sla',
            'parameter': 'Sea level anomaly',
            'timerange': ['0s', '+86399s']
            },
           'Mean_Dynamic_Topo':
           {'productname': 'Mean_Dynamic_Topography',
            'hname': 'mdt',
            'parameter': 'Mean Dynamic Topography',
            'timerange': ['-31d', '+31d']
            },
           'Mean_Sea_Surface':
           {'productname': 'Mean_Sea_Surface',
            'hname': 'Grid_0001',
            'parameter': 'Mean Sea Surface',
            'timerange': ['-31d', '+31d']
            },
           'Sea_Wave_Height':
           {'productname': 'Sea_Wave_Height',
            'hname': 'Grid_0001',
            'parameter': 'Sea Wave Height',
            'timerange': ['0s', '+86399s']
            },
           'Wind':
           {'productname': 'Wind',
            'hname': 'Grid_0001',
            'parameter': 'Wind',
            'timerange': ['0s', '+86399s']
            },
           'Tide':
           {'productname': 'barotropic_tidal_height',
            'hname': 'h',
            'parameter': 'Current tides',
            'timerange': ['-30m', '+30m']
            },
           }


def sea_level_gridded(infile, outdir,
                      vmin=-1., vmax=1.0, vmin_pal=-1., vmax_pal=1.,
                      write_netcdf=False):
    """
    """
    # Read/Process data
    print('Read/Process data')
    ncfile = NCFile(infile)
    if 'id' in ncfile.read_global_attributes():
        l4id = ncfile.read_global_attribute('id')
    elif (re.match(r'^MSL_Map_MERGED_Global_IB_RWT_NoGIA.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Mean_Sea_Level'
        # vmin = -10; vmax = 10; vmin_pal = -10.; vmax_pal = 10
    elif (re.match(r'^dt_global_allsat_madt_h.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Surface_height'
        # vmin = -2; vmax = 2; vmin_pal = -2; vmax_pal = 2
    elif (re.match(r'^dt_global_allsat_msla_h.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Sea_Level_Anomaly'
        # vmin = -0.2; vmax = 0.2; vmin_pal = -0.2; vmax_pal = 0.2
    elif re.match(r'^mdt.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'Mean_Dynamic_Topo'
        # vmin = -1.5; vmax = 1.5; vmin_pal = -1.5; vmax_pal = 1.5
    elif re.match(r'^mss.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'Mean_Sea_Surface'
        # vmin = -80.; vmax = 80.; vmin_pal = -80.; vmax_pal = 80.
    elif (re.match(r'^nrt_merged_mswh.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Sea_Wave_Height'
        # vmin = 0.; vmax = 6.0; vmin_pal = 0.; vmax_pal = 6.0
    elif (re.match(r'^nrt_merged_mwind.*\.nc',
          os.path.basename(infile)) is not None):
        l4id = 'Wind'
        # vmin = 0.; vmax = 20.0; vmin_pal = 0.; vmax_pal = 20.0
    elif re.match(r'^Tide.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'Tide'
        # vmin = -1.5; vmax = 1.5; vmin_pal = -1.5; vmax_pal = 1.50
    else:
        raise Exception('Unknown file.')
    # /TMP
    if l4id in ['Mean_Sea_Level', ]:
        h = ncfile.read_values(L4_MAPS[l4id]['hname'])[::-1, :]
        lon = ncfile.read_values('longitude')[:].astype('float64')
        lat = ncfile.read_values('latitude')[-1:-3:-1].astype('float64')
    elif l4id in ['Mean_Sea_Surface', 'Sea_Wave_Height', 'Wind']:
        h = ncfile.read_values(L4_MAPS[l4id]['hname'])[:, ::-1]
        h = np.transpose(h)
        lon = ncfile.read_values('NbLongitudes')[:].astype('float64')
        lat = ncfile.read_values('NbLatitudes')[-1:-3:-1].astype('float64')
    else:
        h = ncfile.read_values(L4_MAPS[l4id]['hname'])[0, ::-1, :]
        lon = ncfile.read_values('lon')[:].astype('float64')
        lat = ncfile.read_values('lat')[-1:-3:-1].astype('float64')
    lon[lon > 180.] = lon[lon > 180.]-360.
    indsorted = np.argsort(lon)
    lon = lon[indsorted]
    h = h[:, indsorted]
    lon0, dlon, lat0, dlat = lon[0], lon[1]-lon[0], lat[0], lat[1]-lat[0]
    if l4id in ['Mean_Sea_Level', 'Mean_Dynamic_Topo', 'Mean_Sea_Surface']:
        dtime = datetime(2014, 12, 1)
    elif l4id in ['Sea_Wave_Height', 'Wind']:
        dtime = datetime(int(infile[-20:-16]), int(infile[-16:-14]),
                         int(infile[-14:-12]))
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
    if l4id in ['Mean_Sea_Surface', 'Sea_Wave_Height', 'Wind']:
        geolocation['geotransform'] = [lon0-dlon, dlon, 0,
                                       lat0-dlat, 0, dlat]
    elif l4id in ['Tide', ]:
        h = h[:-1, ::]
        geolocation['geotransform'] = [lon0, dlon, 0,
                                       lat0, 0, dlat]
    else:
        geolocation['geotransform'] = [lon0-dlon/2., dlon, 0,
                                       lat0-dlat/2., 0, dlat]
    band = []
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(h, vmin, vmax, out=h)
    array = np.round((h - offset) / scale).astype('uint8')
    array[h.mask] = 255
    colortable = stfmt.format_colortable('matplotlib_jet', vmin=vmin,
                                         vmax=vmax, vmin_pal=vmin_pal,
                                         vmax_pal=vmax_pal)
    band.append({'array': array, 'scale': scale, 'offset': offset,
                 'description': L4_MAPS[l4id]['hname'],
                 'unittype': 'm', 'nodatavalue': 255,
                 'parameter_range': [vmin, vmax], 'colortable': colortable})
    # Write geotiff
    print('Write geotiff')
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
