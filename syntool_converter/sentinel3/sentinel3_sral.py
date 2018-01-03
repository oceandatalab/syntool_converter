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
import numpy
from netCDF4 import Dataset, num2date, date2num
import pyproj
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
import os
import re
import sys
import errno
import tools_for_gcp
from scipy import stats, interpolate

logger = logging.getLogger(__name__)

TIMEFMT = '%Y-%m-%dT%H:%M:%S.%fZ'
cnes_day = datetime(1970,1,1,0,0)
julian_day = datetime(1950,1,1,0,0)

L2_MAPS = {
    'Ku':
    {'productname': 'Sentinel-3_SRAL_20Hz',
     'hname': 'ssha_20_ku',
     'wname': 'swh_ocean_20_ku',
     'sname': 'sig0_ocean_20_ku',
     'lonname': 'lon_20_ku',
     'latname': 'lat_20_ku',
     'timename': 'time_20_ku',
     'parameter': 'S3_SRAL_20Hz',
     'surface': 'surf_type_20_ku'
    },
    '1Hz':
    {'productname': 'Sentinel-3_SRAL_1Hz',
     'hname': 'ssha_01_ku',
     'wname': 'swh_ocean_01_ku',
     'sname': 'sig0_ocean_01_ku',
     'winame': 'wind_speed_alt_01_ku',
     'lonname': 'lon_01',
     'latname': 'lat_01',
     'timename': 'time_01',
     'parameter': 'S3_SRAL_1Hz',
     'surface': 'surf_type_01'
    },
}

def read_orbit(infile, outdir,
               vmin=-0.2, vmax=0.2, vmin_pal=-0.2, vmax_pal=0.2,
               dist_gcp=5, process_var=('ssha', 'swh', 'ws', 'sigma0'),
               keep_empty=False, log_path=None):
    """
    """
    t0 = datetime.utcnow()
    syntool_stats = {}
    if log_path is not None and not os.path.exists(log_path):
        try:
            os.makedirs(log_path)
        except OSError:
            _, e, _ = sys.exc_info()
            if e.errno != errno.EEXIST:
                raise

    listvar=[]
    height = False
    wave = False
    if 'ssha' in process_var:
        height = True
    if 'swh' in process_var:
        wave = True
    if 'ws' in process_var:
        wind = True
    if 'sigma0' in process_var:
        sigma = True
    if height is False and wave is False and wind is False and sigma is False:
       sys.exit('please specify at least on variable to convert')
    L2id = '1Hz' #'Ku'
    # Read variables
    dset = Dataset(os.path.join(infile, 'standard_measurement.nc'), 'r')
    ### For debug purposes
    ind_tmp=0
    lon_0 = dset.variables[L2_MAPS[L2id]['lonname']][ind_tmp:]
    lon_0[lon_0 > 180] = lon_0[lon_0 > 180] - 360
    lat_0 = dset.variables[L2_MAPS[L2id]['latname']][ind_tmp:]
    mask_surface = dset.variables[L2_MAPS[L2id]['surface']][ind_tmp:]
    if height is True:
        ssha_0 = dset.variables[L2_MAPS[L2id]['hname']][ind_tmp:]
        ssha_fill_value = dset.variables[L2_MAPS[L2id]['hname']]._FillValue
        ssha_0[mask_surface > 1] = numpy.nan
    if wave is True:
        swh_0 = dset.variables[L2_MAPS[L2id]['wname']][ind_tmp:]
        swh_fill_value = dset.variables[L2_MAPS[L2id]['wname']]._FillValue
        swh_0[mask_surface > 1] = numpy.nan
    if wind is True:
        ws_0 = dset.variables[L2_MAPS[L2id]['winame']][ind_tmp:]
        ws_fill_value = dset.variables[L2_MAPS[L2id]['winame']]._FillValue
        ws_0[mask_surface > 1] = numpy.nan
    if sigma is True:
        sigma0_0 = dset.variables[L2_MAPS[L2id]['sname']][ind_tmp:]
        sigma0_fill_value = dset.variables[L2_MAPS[L2id]['sname']]._FillValue
        sigma0_0[mask_surface > 1] = numpy.nan

    # var_0 /=  dset.variables[L2_MAPS[L2id]['hname']].scale_factor
    time_0 = dset.variables[L2_MAPS[L2id]['timename']][ind_tmp:]
    time_units = dset.variables[L2_MAPS[L2id]['timename']].units
    dset.close()
    # Trick to deal with continuity in longitude
    dlon = abs(lon_0[1:] - lon_0[:-1])
    lref = lon_0[ numpy.shape(lon_0)[0]/2 ]
    lon_0 = numpy.mod(lon_0 - (lref - 180), 360) + (lref - 180)
    lon_0 = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_0)))

    # Interpolate on land
    dtime = time_0[1:] - time_0[:-1]
    dlon = abs(lon_0[1:] - lon_0[:-1])
    # delta = numpy.median(dtime)
    if len(time_0) > 3:
        delta = stats.mode(dtime)[0][0]
    else:
        logger.warn('orbit too short')
        return

    t_start = datetime.utcnow()
    ndelta = numpy.round(dtime / delta).astype('int')
    if L2id == '1Hz':
        dtime_threshold = 3
    else:
        dtime_threshold = 3 / 20
    ind_dtime = numpy.where(ndelta >= dtime_threshold)[0]
    if len(ind_dtime) == 0:
        time = time_0
        lon = lon_0
        lat = lat_0
        if height is True:
            ssha = ssha_0
        if wave is True:
            swh = swh_0
        if wind is True:
            ws = ws_0
        if sigma is True:
            sigma0 = sigma0_0
    else:
        time = []
        lon = []
        lat = []
        time.append(time_0[:ind_dtime[0] + 1])
        lon.append(lon_0[:ind_dtime[0] + 1])
        lat.append(lat_0[:ind_dtime[0] + 1])
        func_lon = interpolate.interp1d(time_0, lon_0, kind='cubic')
        func_lat = interpolate.interp1d(time_0, lat_0, kind='cubic')
        if height is True:
            ssha = []
            ssha.append(ssha_0[:ind_dtime[0] + 1])
        if wave is True:
            swh = []
            swh.append(swh_0[:ind_dtime[0] + 1])
        if wind is True:
            ws = []
            ws.append(ws_0[:ind_dtime[0] + 1])
        if sigma is True:
            sigma0 = []
            sigma0.append(sigma0_0[:ind_dtime[0] + 1])
        for i in range(len(ind_dtime)):
            time_fill = numpy.linspace(time_0[ind_dtime[i]],
                                       time_0[ind_dtime[i] + 1],
                                       num=((time_0[ind_dtime[i]+1]
                                            - time_0[ind_dtime[i]])/delta +1))
            time.append(time_fill[1:-1])
            lon_tmp = func_lon(time_fill[1:-1])
            lat_tmp = func_lat(time_fill[1:-1])
            lon.append(lon_tmp)
            lat.append(lat_tmp)
            if height is True:
                ssha_fill = numpy.full(numpy.shape(time_fill[1:-1]), numpy.nan)
                ssha.append(ssha_fill)
            if wave is True:
                swh_fill = numpy.full(numpy.shape(time_fill[1:-1]), numpy.nan)
                swh.append(swh_fill)
            if wind is True:
                ws_fill = numpy.full(numpy.shape(time_fill[1:-1]), numpy.nan)
                ws.append(ws_fill)
            if sigma is True:
                sigma0_fill = numpy.full(numpy.shape(time_fill[1:-1]),
                                         numpy.nan)
                sigma0.append(sigma0_fill)
            if i != (len(ind_dtime) - 1):
                slice_ind =  slice(ind_dtime[i] + 1, ind_dtime[i + 1] + 1)
                time.append(time_0[slice_ind])
                lon.append(lon_0[slice_ind])
                lat.append(lat_0[slice_ind])
                if height is True:
                    ssha.append(ssha_0[slice_ind])
                if wave is True:
                    swh.append(swh_0[slice_ind])
                if wind is True:
                    ws.append(ws_0[slice_ind])
                if sigma is True:
                    sigma0.append(sigma0_0[slice_ind])
            else:
                time.append(time_0[ind_dtime[i] + 1:])
                lon.append(lon_0[ind_dtime[i] + 1:])
                lat.append(lat_0[ind_dtime[i] + 1:])
                if height is True:
                    ssha.append(ssha_0[ind_dtime[i] + 1:])
                if wave is True:
                    swh.append(swh_0[ind_dtime[i] + 1:])
                if wind is True:
                    ws.append(ws_0[ind_dtime[i] + 1:])
                if sigma is True:
                    sigma0.append(sigma0_0[ind_dtime[i] + 1:])
        time = numpy.concatenate(time, axis=0)
        lon = numpy.concatenate(lon, axis=0)
        lat = numpy.concatenate(lat, axis=0)
        if height is True:
            ssha = numpy.concatenate(ssha, axis=0)
            ssha[ssha == ssha_fill_value] = numpy.nan
        if wave is True:
            swh = numpy.concatenate(swh, axis=0)
            swh[swh == swh_fill_value] = numpy.nan
        if wind is True:
            ws = numpy.concatenate(ws, axis=0)
            ws[ws == ws_fill_value] = numpy.nan
        if sigma is True:
            sigma0 = numpy.concatenate(sigma0, axis=0)
            sigma0[sigma0 == sigma0_fill_value] = numpy.nan
    t_stop = datetime.utcnow()
    syntool_stats['gapfill_computation'] = (t_stop - t_start).total_seconds()
    lon = lon - numpy.floor((numpy.min(lon) + 180.) / 360.) * 360.
    ntime = numpy.shape(time)[0]
    if height is True:
        nan_mask = numpy.isnan(ssha)
        nan_ind = numpy.where(nan_mask)
        ssha[nan_ind] = 0
        mask_gap_ssha = (nan_mask | (abs(ssha) > 50))
        ssha[nan_ind] = numpy.nan  # Restore nan values
        if not mask_gap_ssha.all():
            _min = numpy.nanmin(ssha[numpy.where(~mask_gap_ssha)])
            _max = numpy.nanmax(ssha[numpy.where(~mask_gap_ssha)])
            syntool_stats['ssha'] = {'min': _min, 'max': _max}
        listvar.append({'var': ssha, 'mask_gap': mask_gap_ssha,
                        'range': [-0.4, 0.4], 'parameter': 'SSHA',
                        'palette':'matplotlib_jet',
                        'range_pal': [-0.4, 0.4]})
    if wave is True:
        nan_mask = numpy.isnan(swh)
        nan_ind = numpy.where(nan_mask)
        swh[nan_ind] = 0
        mask_gap_swh = (nan_mask | (abs(swh) > 50))
        swh[nan_ind] = numpy.nan  # Restore nan values
        if not mask_gap_swh.all():
            _min = numpy.nanmin(swh[numpy.where(~mask_gap_swh)])
            _max = numpy.nanmax(swh[numpy.where(~mask_gap_swh)])
            syntool_stats['swh'] = {'min': _min, 'max': _max}
        listvar.append({'var': swh, 'mask_gap': mask_gap_swh,
                        'range': [0., 8.], 'parameter': 'SWH',
                        'palette':'matplotlib_jet',
                        'range_pal': [0., 8.]})
    if wind is True:
        nan_mask = numpy.isnan(ws)
        nan_ind = numpy.where(nan_mask)
        ws[nan_ind] = 0
        mask_gap_ws = (nan_mask | (abs(ws) > 200))
        ws[nan_ind] = numpy.nan  # Restore nan values
        if not mask_gap_ws.all():
            _min = numpy.nanmin(ws[numpy.where(~mask_gap_ws)])
            _max = numpy.nanmax(ws[numpy.where(~mask_gap_ws)])
            syntool_stats['ws'] = {'min': _min, 'max': _max}
        listvar.append({'var': ws, 'mask_gap': mask_gap_ws,
                        'range': [0., 25.], 'parameter': 'Ws',
                        'palette':'noaa_wind',
                        'range_pal': [0., 25.]})
    if sigma is True:
        nan_mask = numpy.isnan(sigma0)
        nan_ind = numpy.where(nan_mask)
        sigma0[nan_ind] = 0
        mask_gap_sigma0 = (nan_mask | (sigma0 > 100))
        sigma0[nan_ind] = numpy.nan  # Restore nan values
        if not mask_gap_sigma0.all():
            _min = numpy.nanmin(sigma0[numpy.where(~mask_gap_sigma0)])
            _max = numpy.nanmax(sigma0[numpy.where(~mask_gap_sigma0)])
            syntool_stats['sigma0'] = {'min': _min, 'max': _max}
        listvar.append({'var': sigma0, 'mask_gap': mask_gap_sigma0,
                        'range': [5., 25.], 'parameter': 'Sigma0',
                        'palette':'matplotlib_jet',
                        'range_pal': [5., 25.]})
    start_time = num2date(time[0], time_units)
    end_time = num2date(time[-1], time_units)
    # NOTE : lon/lat must be continuous even if crossing dateline
    # (ie. no [-180,180] clipping)
    # Make GCPs (mimic a swath of arbitrary width in lon/lat, here ~5km)
    # gcps = tools_for_gcp.make_gcps_v1(lon, lat, dist_gcp=dist_gcp)
    gcps = tools_for_gcp.make_gcps_v2(lon, lat,
                                      dist_gcp=dist_gcp)
    file_name, _ = os.path.splitext(os.path.basename(os.path.normpath(infile)))
    for i in range(len(listvar)):
        metadata = {}
        (dtime, time_range) = stfmt.format_time_and_range(start_time,
                                                          end_time,
                                                          units='s')
        metadata['product_name'] = (L2_MAPS[L2id]['productname']
                                    + '_' + listvar[i]['parameter'])
        metadata['name'] = file_name
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
        #ssha[numpy.where(~(numpy.isfinite(mask_gap)))] = 255
        var = listvar[i]['var']
        mask_gap = listvar[i]['mask_gap']
        var[numpy.where(mask_gap)] = 255
        array = numpy.clip(numpy.round((var - offset) / scale),
                           0, 254).astype('uint8')
        array[numpy.where(mask_gap)] = 255
        array = array[:, numpy.newaxis]
        if not keep_empty and numpy.all(array == 255):
            parameter = listvar[i]['parameter']
            logger.warn('No valid values in this dataset for ' \
                        '{}'.format(parameter))
            logger.warn('Skipped.')
            continue
        colortable = stfmt.format_colortable(listvar[i]['palette'],
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

    syntool_stats['total_time'] = (datetime.utcnow() - t0).total_seconds()
    if log_path is not None:
        import json
        stats_path = os.path.join(log_path, '{}.json'.format(file_name))
        with open(stats_path, 'w') as f:
            json.dump(syntool_stats, f)
