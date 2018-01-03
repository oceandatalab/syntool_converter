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

import logging
import numpy as np
from netCDF4 import Dataset, num2date, date2num
import pyproj
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
import os
import re
import sys
import tools_for_gcp
from scipy import stats, interpolate

logger = logging.getLogger(__name__)

TIMEFMT = '%Y-%m-%dT%H:%M:%S.%fZ'
cnes_day = datetime(1970,1,1,0,0)
julian_day = datetime(1950,1,1,0,0)

L2_MAPS = {
    # X_TRACK GoM SLA
    'X_TRACK_GOM_SLA':
    {'productname': 'X_TRACK_GOM_SLA',
     'hname': 'SLA',
     'lonname': 'lon',
     'latname': 'lat',
     'parameter': 'X_TRACK SLA',
     },
    # X_TRACK GoM SLA
    'X_TRACK_MED_SLA':
    {'productname': 'X_TRACK_MED_SLA',
     'hname': 'SLA',
     'lonname': 'lon',
     'latname': 'lat',
     'parameter': 'X_TRACK SLA',
    },
    # X_TRACK GoM SLA
    'X_TRACK_NEA_SLA':
    {'productname': 'X_TRACK_NEA_SLA',
     'hname': 'SLA',
     'lonname': 'lon',
     'latname': 'lat',
     'parameter': 'X_TRACK SLA',
    },
        }


def read_orbit(infile, outdir,
               vmin=-0.2, vmax=0.2, vmin_pal=-0.2, vmax_pal=0.2,
               dist_gcp=25.5, keep_empty=False):
    """
    """
    if (re.match('.*SLAXT_filt_NorthEastAtlantic*.', os.path.basename(infile))
          is not None):
        L2id = 'X_TRACK_NEA_SLA'
    elif (re.match('.*SLAXT_filt_GoMCaribbean*.', os.path.basename(infile))
          is not None):
        L2id = 'X_TRACK_GOM_SLA'
    elif (re.match('.*SLAXT_filt_MediterraneanSea*.', os.path.basename(infile))
          is not None):
        L2id = 'X_TRACK_MED_SLA'
    else: 
        logger.warn('Unknown id for file {}'.format(os.path.basename(infile)))
        # vmin = -0.2 ; vmax = 0.2 ; vmin_pal = -0.2 ; vmax_pal = 0.2

    # Read variables
    dset = Dataset(infile)
    ### For debug purposes
    ind_tmp=0
    lon_0 = dset.variables[L2_MAPS[L2id]['lonname']][ind_tmp:]
    lon_0[lon_0 > 180] = lon_0[lon_0 > 180] - 360
    lat_0 = dset.variables[L2_MAPS[L2id]['latname']][ind_tmp:]
    ssha_0 = dset.variables[L2_MAPS[L2id]['hname']][ind_tmp:]
    ssha_fill_value = dset.variables[L2_MAPS[L2id]['hname']]._FillValue

    # var_0 /=  dset.variables[L2_MAPS[L2id]['hname']].scale_factor
    time_0 = dset.variables['time'][ind_tmp:]
    time_units = dset.variables['time'].units
    dset.close()
    # Trick to deal with continuity in longitude
    dlon = abs(lon_0[1:] - lon_0[:-1])
    lref = lon_0[ np.shape(lon_0)[0]/2 ]
    lon_0 = np.mod(lon_0 - (lref - 180), 360) + (lref - 180)
    lon_0 = np.rad2deg(np.unwrap(np.deg2rad(lon_0)))

    # Interpolate on land
    dtime = time_0[1:] - time_0[:-1]
    dlon = abs(lon_0[1:] - lon_0[:-1])
    # delta = np.median(dtime)
    if len(time_0) > 3:
        delta = stats.mode(dtime)[0][0]
    else:
        sys.exit('orbit too short')
    ndelta = np.round(dtime / delta).astype('int')
    ind_dtime = np.where(ndelta >= 5)[0]
    if len(ind_dtime) == 0:
        time = time_0
        ssha = ssha_0
        lon = lon_0
        lat = lat_0
    else:
        time = time_0[:ind_dtime[0] + 1]
        ssha = ssha_0[:ind_dtime[0] + 1]
        for i in range(len(ind_dtime)):
           # time_fill = np.linspace(time_0[ind_dtime[i]], time_0[ind_dtime[i] + 1],
           #                         num=(time_0[ind_dtime[i]+1]
           #                         - time_0[ind_dtime[i]])/delta)
            time_fill = np.linspace(time_0[ind_dtime[i]],
                                    time_0[ind_dtime[i] + 1],
                                    num=ndelta[ind_dtime[i]], endpoint=False)
            ssha_fill = np.zeros(np.shape(time_fill[1:])) * np.nan
            ssha = np.hstack([ssha, ssha_fill])
            time = np.hstack([time, time_fill[1:]])
            if i != (len(ind_dtime) - 1):
                time = np.hstack([time, time_0[ind_dtime[i] + 1:ind_dtime[i + 1]
                                               + 1]])
                ssha = np.hstack([ssha, ssha_0[ind_dtime[i]
                                        + 1:ind_dtime[i + 1] + 1]])
            else:
                time = np.hstack([time, time_0[ind_dtime[i] + 1:]])
                ssha = np.hstack([ssha, ssha_0[ind_dtime[i] + 1:]])
        func = interpolate.interp1d(time_0, lon_0, kind='linear')
        lon = func(time)
        func = interpolate.interp1d(time_0, lat_0, kind='linear')
        lat = func(time)
        ssha[ssha == ssha_fill_value] = np.nan

    lon = lon - np.floor((np.min(lon) + 180.) / 360.) * 360.
    ntime = np.shape(time)[0]
    mask_gap_ssha = np.isnan(ssha)
    start_time = num2date(time[0], time_units)
    end_time = num2date(time[-1], time_units)
    # NOTE : lon/lat must be continuous even if crossing dateline
    # (ie. no [-180,180] clipping)
    # Make GCPs (mimic a swath of arbitrary width in lon/lat, here ~5km)
    # gcps = tools_for_gcp.make_gcps_v1(lon, lat, dist_gcp=dist_gcp)
    gcps = tools_for_gcp.make_gcps_v2(lon, lat,
                                      dist_gcp=dist_gcp)
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time,
                                                      end_time,
                                                      units='s')
    metadata['product_name'] = (L2_MAPS[L2id]['productname'])
    metadata['name'] = (os.path.splitext(os.path.basename(infile))[0])
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['gcps'] = stfmt.format_gdalgcps(*gcps)
    band = []
    #vmin = listvar[i]['range'][0]
    #vmax = listvar[i]['range'][1]
    #vmin_pal = listvar[i]['range_pal'][0]
    #vmax_pal = listvar[i]['range_pal'][1]
    scale = (vmax - vmin) / 254.
    offset = vmin
        # Avoid warnings caused by NaN values
        #ssha[np.where(~(np.isfinite(mask_gap)))] = 255
    var = ssha
    var[np.where(mask_gap_ssha)] = 255
    array = np.clip(np.round((var - offset) / scale), 0, 254).astype('uint8')
    array[np.where(mask_gap_ssha)] = 255
    array = array[:, np.newaxis]
    if not keep_empty and np.all(array == 255):
        parameter = L2_MAPS[L2id]['parameter']
        logger.warn('No valid values in this dataset for ' \
                    '{}'.format(parameter))
        logger.warn('Skipped.')
    colortable = stfmt.format_colortable('matplotlib_jet',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array': array, 'scale': scale, 'offset': offset,
                 'description':  L2_MAPS[L2id]['parameter'], 'unittype': 'm',
                 'nodatavalue': 255, 'parameter_range': [vmin, vmax],
                 'colortable': colortable})
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    # Check GCPs
    # tools_for_gcp.check_gcps(tifffile, lon, lat, gcps[0], gcps[1])
