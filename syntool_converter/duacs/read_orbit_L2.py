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
    # Jason 2 GDR
    'J2_GDR':
    {'productname': 'Jason_2_GDR',
     'hname': 'ssha',
     'wname': 'swh_ku',
     'lonname': 'lon',
     'latname': 'lat',
     'parameter': 'Jason 2 GDR',
     'surface': 'surface_type'
    },
    # Saral GDR
    'SRL_GDR':
    {'productname': 'Saral_GDR',
     'hname': 'ssha',
     'wname': 'swh',
     'lonname': 'lon',
     'latname': 'lat',
     'parameter': 'Saral GDR',
     'surface': 'surface_type'
    },
        }


def read_orbit(infile, outdir,
               vmin=-0.2, vmax=0.2, vmin_pal=-0.2, vmax_pal=0.2,
               dist_gcp=12.5, process_var=['ssha', 'swh'], keep_empty=False):
    """
    """
    listvar=[]
    height = False
    wave = False
    if 'ssha' in process_var:
        height = True
    if 'swh' in process_var:
        wave = True
    if height is False and wave is False:
       sys.exit('please specify at least on variable to convert')
    if re.match(r'^SRL_.*\.nc', os.path.basename(infile)) is not None:
        L2id = 'SRL_GDR'
    elif re.match(r'^JA2_.*\.nc', os.path.basename(infile)) is not None:
        L2id = 'J2_GDR'

    # Read variables
    dset = Dataset(infile, 'r')
    ### For debug purposes
    ind_tmp=0
    lon_0 = dset.variables[L2_MAPS[L2id]['lonname']][ind_tmp:]
    lon_0[lon_0 > 180] = lon_0[lon_0 > 180] - 360
    lat_0 = dset.variables[L2_MAPS[L2id]['latname']][ind_tmp:]
    mask_surface = dset.variables[L2_MAPS[L2id]['surface']][ind_tmp:]
    if height is True:
        ssha_0 = dset.variables[L2_MAPS[L2id]['hname']][ind_tmp:]
        ssha_fill_value = dset.variables[L2_MAPS[L2id]['hname']]._FillValue
        ssha_0[mask_surface > 1] = np.nan
    if wave is True:
        swh_0 = dset.variables[L2_MAPS[L2id]['wname']][ind_tmp:]
        swh_fill_value = dset.variables[L2_MAPS[L2id]['wname']]._FillValue
        swh_0[mask_surface > 1] = np.nan

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
    ind_dtime = np.where(ndelta >= 2)[0]
    if len(ind_dtime) == 0:
        time = time_0
        ssha = ssha_0
        lon = lon_0
        lat = lat_0
        if wave is True:
            swh = swh_0
    else:
        time = time_0[:ind_dtime[0] + 1]
        if height is True:
            ssha = ssha_0[:ind_dtime[0] + 1]
        if wave is True:
            swh = swh_0[:ind_dtime[0] + 1]
        for i in range(len(ind_dtime)):
            time_fill = np.linspace(time_0[ind_dtime[i]], time_0[ind_dtime[i] + 1],
                                    num=(time_0[ind_dtime[i]+1]
                                    - time_0[ind_dtime[i]])/delta)
            if height is True:
                ssha_fill = np.zeros(np.shape(time_fill[1:])) * np.nan
                ssha = np.hstack([ssha, ssha_fill])
            if wave is True:
                swh_fill = np.zeros(np.shape(time_fill[1:])) * np.nan
                swh = np.hstack([swh, swh_fill])
            time = np.hstack([time, time_fill[1:]])
            if i != (len(ind_dtime) - 1):
                time = np.hstack([time, time_0[ind_dtime[i] + 1:ind_dtime[i + 1]
                                               + 1]])
                if height is True:
                    ssha = np.hstack([ssha, ssha_0[ind_dtime[i]
                                            + 1:ind_dtime[i + 1] + 1]])
                if wave is True:
                    swh = np.hstack([swh, swh_0[ind_dtime[i]
                                            + 1:ind_dtime[i + 1] + 1]])
            else:
                time = np.hstack([time, time_0[ind_dtime[i] + 1:]])
                if height is True:
                    ssha = np.hstack([ssha, ssha_0[ind_dtime[i] + 1:]])
                if wave is True:
                    swh = np.hstack([swh, swh_0[ind_dtime[i] + 1:]])
        func = interpolate.interp1d(time_0, lon_0, kind='quadratic')
        lon = func(time)
        func = interpolate.interp1d(time_0, lat_0, kind='quadratic')
        lat = func(time)
        if height is True:
            ssha[ssha == ssha_fill_value] = np.nan
        if wave is True:
            swh[swh == swh_fill_value] = np.nan
    lon = lon - np.floor((np.min(lon) + 180.) / 360.) * 360.
    ntime = np.shape(time)[0]
    if height is True:
        mask_gap_ssha = np.isnan(ssha)
        listvar.append({'var': ssha, 'mask_gap': mask_gap_ssha,
                    'range': [-0.4, 0.4], 'parameter': 'SSHA',
                    'range_pal': [-0.4, 0.4]})
    if wave is True:
        mask_gap_swh = np.isnan(swh)
        listvar.append({'var': swh, 'mask_gap': mask_gap_swh,
                        'range': [0., 8.], 'parameter': 'SWH',
                        'range_pal': [0., 8.]})
    start_time = num2date(time[0], time_units)
    end_time = num2date(time[-1], time_units)
    # NOTE : lon/lat must be continuous even if crossing dateline
    # (ie. no [-180,180] clipping)
    # Make GCPs (mimic a swath of arbitrary width in lon/lat, here ~5km)
    # gcps = tools_for_gcp.make_gcps_v1(lon, lat, dist_gcp=dist_gcp)
    gcps = tools_for_gcp.make_gcps_v2(lon, lat,
                                      dist_gcp=dist_gcp)
    for i in range(len(listvar)):
        metadata = {}
        (dtime, time_range) = stfmt.format_time_and_range(start_time,
                                                          end_time,
                                                          units='s')
        metadata['product_name'] = (L2_MAPS[L2id]['productname']
                                    + '_' + listvar[i]['parameter'])
        metadata['name'] = (os.path.splitext(os.path.basename(infile))[0])
        metadata['datetime'] = dtime
        metadata['time_range'] = time_range
        geolocation = {}
        geolocation['projection'] = stfmt.format_gdalprojection()
        geolocation['gcps'] = stfmt.format_gdalgcps(*gcps)
        band = []
        vmin = listvar[i]['range'][0]
        vmax = listvar[i]['range'][1]
        vmin_pal = listvar[i]['range_pal'][0]
        vmax_pal = listvar[i]['range_pal'][1]
        scale = (vmax - vmin) / 254.
        offset = vmin
        # Avoid warnings caused by NaN values
        #ssha[np.where(~(np.isfinite(mask_gap)))] = 255
        var = listvar[i]['var']
        mask_gap = listvar[i]['mask_gap']
        var[np.where(mask_gap)] = 255
        array = np.clip(np.round((var - offset) / scale), 0, 254).astype('uint8')
        array[np.where(mask_gap)] = 255
        array = array[:, np.newaxis]
        if not keep_empty and np.all(array == 255):
            parameter = listvar[i]['parameter']
            logger.warn('No valid values in this dataset for ' \
                        '{}'.format(parameter))
            logger.warn('Skipped.')
            continue
        colortable = stfmt.format_colortable('matplotlib_jet',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
        band.append({'array': array, 'scale': scale, 'offset': offset,
                     'description': listvar[i]['parameter'], 'unittype': 'm',
                     'nodatavalue': 255, 'parameter_range': [vmin, vmax],
                     'colortable': colortable})
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    # Check GCPs
    # tools_for_gcp.check_gcps(tifffile, lon, lat, gcps[0], gcps[1])
