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
    'CLS-L4-CUReul_hs-ALT_SUM-v01.0':
    {
        'productname': 'GlobCurrent_L4_total_hs',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CUReul_15m-ALT_SUM-v01.0':
    {
        'productname': 'GlobCurrent_L4_total_15m',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CURgeo_0m-ALT_OI-v01.0':
    {
        'productname': 'GlobCurrent_L4_geostrophic',
        'uname': 'eastward_geostrophic_current_velocity',
        'vname': 'northward_geostrophic_current_velocity',
        'timerange': ['-12h', '+12h']
    },
    'CLS-L4-CURekm_hs-ERAWS_EEM-v01.0':
    {
        'productname': 'GlobCurrent_L4_ekman_hs',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CURekm_15m-ERAWS_EEM-v01.0':
    {
        'productname': 'GlobCurrent_L4_ekman_15m',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    # New format v02
    'CLS-L4-CUReul_hs-ALT_SUM-v02.0':
    {
        'productname': 'GlobCurrent_L4_total_hs',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CUReul_15m-ALT_SUM-v02.0':
    {
        'productname': 'GlobCurrent_L4_total_15m',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CURgeo_0m-ALT_OI-v02.0':
    {
        'productname': 'GlobCurrent_L4_geostrophic',
        'uname': 'eastward_geostrophic_current_velocity',
        'vname': 'northward_geostrophic_current_velocity',
        'timerange': ['-12h', '+12h']
    },
    'CLS-L4-CURekm_hs-ERAWS_EEM-v02.0':
    {
        'productname': 'GlobCurrent_L4_ekman_hs',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CURekm_15m-ERAWS_EEM-v02.0':
    {
        'productname': 'GlobCurrent_L4_ekman_15m',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    # New format v03
    'CLS-L4-CUReul_hs-ALT_SUM-v03.0':
    {
        'productname': 'GlobCurrent_L4_total_hs',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CUReul_15m-ALT_SUM-v03.0':
    {
        'productname': 'GlobCurrent_L4_total_15m',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CURgeo_0m-ALT_OI-v03.0':
    {
        'productname': 'GlobCurrent_L4_geostrophic',
        'uname': 'eastward_geostrophic_current_velocity',
        'vname': 'northward_geostrophic_current_velocity',
        'timerange': ['-12h', '+12h']
    },
    'CLS-L4-CURekm_hs-ERAWS_EEM-v03.0':
    {
        'productname': 'GlobCurrent_L4_ekman_hs',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'CLS-L4-CURekm_15m-ERAWS_EEM-v03.0':
    {
        'productname': 'GlobCurrent_L4_ekman_15m',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-90m', '+90m']
    },
    'IFREMER-L4-CURstk_hs-WW3_IFR-v01.0':
    {
        'productname': 'GlobCurrent_L4_stokes_drift',
        'uname': 'eastward_stokes_drift_velocity',
        'vname': 'northward_stokes_drift_velocity',
        'timerange': ['-90m', '+90m']
    },

    # Mediterranean product
    'CLS-L4-CURekm_15m-NOAAW_MED_EEM-v03.0':
    {
        'productname': 'GlobCurrent_L4_ekman_15m_med',
        'uname': 'eastward_ekman_current_velocity',
        'vname': 'northward_ekman_current_velocity',
        'timerange': ['-180m', '+180m']
    },
    'CLS-L4-CUReul_15m-ALT_MED_SUM-v03.0':
    {
        'productname': 'GlobCurrent_L4_total_15m_med',
        'uname': 'eastward_eulerian_current_velocity',
        'vname': 'northward_eulerian_current_velocity',
        'timerange': ['-180m', '+180m']
    },
    'CLS-L4-CURgeo_0m-ALT_MED_OI-v03.0':
    {
        'productname': 'GlobCurrent_L4_geostrophic_med',
        'uname': 'eastward_geostrophic_current_velocity',
        'vname': 'northward_geostrophic_current_velocity',
        'timerange': ['-12h', '+12h']
    },

    #NRT product
    'CLS-L4-CURgeo_0m-ALT_OI_NRT-v03.0':
    {
       'productname': 'GlobCurrent_L4_geostrophic_nrt',
       'uname': 'u',
       'vname': 'v',
       'timerange': ['-12h', '+12h']
    },
    'CLS-L4-CURekm_15m-ERAWS_EEM_NRT-v03.0':
    {
       'productname': 'GlobCurrent_L4_ekman_15m_nrt',
       'uname': 'u',
       'vname': 'v',
       'timerange': ['-12h', '+12h']
    },
    'CLS-L4-CURekm_15m-ERAWS_EEM_NRT-6h-v03.0':
    {
       'productname': 'GlobCurrent_L4_ekman_15m_nrt',
       'uname': 'u',
       'vname': 'v',
       'timerange': ['-180m', '+180m']
    },
   'CLS-L4-CUReul_15m-ALT_SUM_NRT-v03.0':
    {
       'productname': 'GlobCurrent_L4_total_15m_nrt',
       'uname': 'u',
       'vname': 'v',
       'timerange': ['-12h', '+12h']
    },
    # Old format
    'GC_L4_CUR_GLO_010_TOTHS':
    {
        'productname': 'GlobCurrent_L4_total_hs',
        'uname': 'eastward_total_current_velocity_hs',
        'vname': 'northward_total_current_velocity_hs',
        'timerange': ['-90m', '+90m']
    },
    'GC_L4_CUR_GLO_010_TOT15':
    {
        'productname': 'GlobCurrent_L4_total_15m',
        'uname': 'eastward_total_current_velocity_15m',
        'vname': 'northward_total_current_velocity_15m',
        'timerange': ['-90m', '+90m']
    },
    'GC_L4_CUR_GLO_010_GEOST':
    {
        'productname': 'GlobCurrent_L4_geostrophic',
        'uname': 'eastward_geostrophic_velocity',
        'vname': 'northward_geostrophic_velocity',
        'timerange': ['-12h', '+12h']
    },
    'GC_L4_CUR_GLO_010_EKHS':
    {
        'productname': 'GlobCurrent_L4_ekman_hs',
        'uname': 'eastward_ekman_velocity_hs',
        'vname': 'northward_ekman_velocity_hs',
        'timerange': ['-90m', '+90m']
    },
    'GC_L4_CUR_GLO_010_EK15':
    {
        'productname': 'GlobCurrent_L4_ekman 15m',
        'uname': 'eastward_ekman_velocity_15m',
        'vname': 'northward_ekman_velocity_15m',
        'timerange': ['-90m', '+90m']
    },
    'GC_MOD_STK_GLO_010_WW3':
    {
        'productname': 'GlobCurrent_L4_stokes_drift',
        'uname': 'eastward_stokes_drift',
        'vname': 'northward_stokes_drift',
        'timerange': ['-90m', '+90m']
    },
    'GC_MOD_TIDE_GLO_010_FES2012':
    {
        'productname': 'GlobCurrent_L4_tidal',
        'uname': 'eastward_tidal_current',
        'vname': 'northward_tidal_current',
        'timerange': ['-90m', '+90m']
    },
    # For testing GlobCurrent reprocessing
    '':
    {
        'productname': 'GlobCurrent_L4_intermediate',
        'uname': 'interpolated_u',
        'vname': 'interpolated_v',
        'timerange': ['-90m', '+90m'] # may be wrong
    },
    # For CLS surcouf NRT (geostrophic)
    'CourantGeostr':
    {
        'productname': 'GlobCurrent_L4_geostrophic',
        'uname': 'u',
        'vname': 'v',
        'timerange': ['-12h', '+12h']
    }
}


# def tmp():
#     import glob
#     patt = '/local/home/data/ancillary/globcurrent/global_010_deg/*/2012/245/*.nc'
#     paths = glob.glob(patt)
#     for path in paths:
#         print '\n'+os.path.basename(path)
#         ncfile = NCFile(path)
#         print ncfile.read_global_attribute('id')
#         print ncfile.read_global_attribute('title')
#         print ncfile.get_fieldnames()


def globcurrent_l4(infile, outdir,
                   vmin=0., vmax=5.08, vmin_pal=0., vmax_pal=2.,
                   write_netcdf=False):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    ncfile = NCFile(infile)
    if 'id' in ncfile.read_global_attributes():
        l4id = ncfile.read_global_attribute('id')
    elif re.match(r'^CourantGeostr_.*\.nc', os.path.basename(infile)) is not None:
        l4id = 'CourantGeostr'
    else:
        raise Exception('Unknown GlobCurrent L4 file.')
    # TMP : bug in GlobCurrent products v02.0
    # Correct bug in first release of v02.0
    # if l4id == 'L4-CURekm_15m-ERAWS_EEM-v02.0' or \
    #    l4id == 'L4-CURekm_hs-ERAWS_EEM-v02.0':
    #     if 'L4-CURekm_hs-ERAWS_EEM-v02.0' in os.path.basename(infile):
    #         l4id = 'CLS-L4-CURekm_hs-ERAWS_EEM-v02.0'
    #     elif 'L4-CURekm_15m-ERAWS_EEM-v02.0' in os.path.basename(infile):
    #         l4id = 'CLS-L4-CURekm_15m-ERAWS_EEM-v02.0'
    #     else:
    #         raise Exception()
    # Check bug is not anymore in second release of v02.0
    if 'v03.0' in l4id:
        if l4id not in L4_MAPS:
            raise Exception('Unknown l4id.')
        if l4id[4:] not in os.path.basename(infile):
            raise Exception('Inconsistency between filename and l4id.')
        if 'CLS-L4-CURekm_15m-ERAWS_EEM_NRT-v03.0' == l4id and re.match('\d{14}-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM_NRT-v03.0-fv01.0.nc', os.path.basename(infile)) is not None:
            l4id = 'CLS-L4-CURekm_15m-ERAWS_EEM_NRT-6h-v03.0'
    # /TMP
    ucur = ncfile.read_values(L4_MAPS[l4id]['uname'])[0, ::-1, :]
    vcur = ncfile.read_values(L4_MAPS[l4id]['vname'])[0, ::-1, :]
    ufattr = ncfile.read_field_attributes(L4_MAPS[l4id]['uname'])
    vfattr = ncfile.read_field_attributes(L4_MAPS[l4id]['vname'])
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
    dtime_units = ncfile.read_field('time').units
    dtime = num2date(ncfile.read_values('time')[0], dtime_units)
    # rundtime = ncfile.read_global_attribute('date_modified')
    # rundtime = datetime.strptime(rundtime, '%Y%m%dT%H%M%SZ')
    ncfile.close()
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = L4_MAPS[l4id]['productname']
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = stfmt.format_time(dtime)
    #metadata['time_range'] = ['-90m', '+90m']
    metadata['time_range'] = L4_MAPS[l4id]['timerange']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'GlobCurrent'
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
                     'description':'current u', 'unittype':'m s-1',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax]})
        if l4id == 'CourantGeostr':
            band[0]['name'] = L4_MAPS['CLS-L4-CURgeo_0m-ALT_OI-v02.0']['uname']
        else:
            band[0]['name'] = L4_MAPS[l4id]['uname']
        if 'long_name' in ufattr:
            band[0]['long_name'] = ufattr['long_name']
        v = np.clip(vcur.data, vmin, vmax)
        array = np.round((v - offset) / scale).astype('uint8')
        array[mask] = 255
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'current v', 'unittype':'m s-1',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax]})
        if l4id == 'CourantGeostr':
            band[1]['name'] = L4_MAPS['CLS-L4-CURgeo_0m-ALT_OI-v02.0']['vname']
        else:
            band[1]['name'] = L4_MAPS[l4id]['vname']
        if 'long_name' in vfattr:
            band[1]['long_name'] = vfattr['long_name']
        # Write
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        metadata['spatial_resolution'] = min([abs(dlat), abs(dlon)]) * 111000.
        dgcps = np.round(1. / np.abs(np.array([dlat, dlon]))).astype('int')
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'grid_lonlat',
                           dgcps=dgcps)
