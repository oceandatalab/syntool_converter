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
from syntool_converter.utils.write_idf1d import write_netcdf_1d
logger = logging.getLogger(__name__)

TIMEFMT = '%Y-%m-%dT%H:%M:%S.%fZ'
cnes_day = datetime(1970,1,1,0,0)
julian_day = datetime(1950,1,1,0,0)

L3_MAPS = {
          # Jason ADT
          'J2_ADT':
          {'productname': 'Jason_2_ADT',
           'hname': 'ADT',
           'parameter': 'Jason 2 ADT',
           },
          # SARAL ADT
          'AL_ADT':
          {'productname': 'Saral_ADT',
           'hname': 'ADT',
           'parameter': 'Saral ADT',
           },
          # Jason SLA
          'J2_SLA':
          {'productname': 'Jason_2_SLA',
           'hname': 'SLA',
           'parameter': 'Jason 2 SLA',
           },
          # SARAL SLA
          'AL_SLA':
          {'productname': 'Saral_SLA',
           'hname': 'SLA',
           'parameter': 'Saral SLA',
           },
          }


def read_orbit(infile, outdir,
               vmin=-0.2, vmax=0.2, vmin_pal=-0.2, vmax_pal=0.2,
               dist_gcp=None, write_netcdf=False):
    """
    """
    if (re.match(r'^dt_global_j2_adt_vfec_.*\.nc', os.path.basename(infile))
          is not None):
        L3id = 'J2_ADT'
        # vmin = -0.2; vmax = 0.2 ; vmin_pal = -0.2 ; vmax_pal = 0.2
    elif (re.match(r'^dt_global_al_adt_vfec.*\.nc', os.path.basename(infile))
          is not None):
        L3id = 'AL_ADT'
        # vmin = -0.2 ; vmax = 0.2 ; vmin_pal = -0.2 ; vmax_pal = 0.2
    elif (re.match(r'^dt_global_j2_sla_vfec.*\.nc', os.path.basename(infile))
          is not None):
        L3id = 'J2_SLA'
        # vmin = -0.2 ; vmax = 0.2 ; vmin_pal = -0.2 ; vmax_pal = 0.2
    elif (re.match(r'^dt_global_al_sla_vfec.*\.nc', os.path.basename(infile))
          is not None):
        L3id = 'AL_SLA'
        # vmin = -0.2 ; vmax = 0.2 ; vmin_pal = -0.2 ; vmax_pal = 0.2

    # Read
    dset = Dataset(infile)
    lon_all = dset.variables['longitude'][:]
    lat_all = dset.variables['latitude'][:]
    var_all = dset.variables[L3_MAPS[L3id]['hname']][:]
    time_all = dset.variables['time'][:]
    ipass_all = dset.variables['track'][:]
    cycle_all = dset.variables['cycle'][:]
    time_units = dset.variables['time'].units
    var_fill_value = dset.variables[L3_MAPS[L3id]['hname']]._FillValue
    var_all[var_all == var_fill_value] = np.nan
    dset.close()

    # Detect passes index
    dipass = ipass_all[1:] - ipass_all[:-1]
    ind_ipass = np.where(dipass != 0)[0]
    ind_ipass += 1
    ind_ipass = np.hstack([0, ind_ipass, len(ipass_all)])

    # Loop on all passes and compute geotiff for each pass
    for i in range(len(ind_ipass)-1):
        lon_0 = lon_all[ind_ipass[i]:ind_ipass[i + 1]]
        # lon_0 = np.mod(lon_0 + 180, 360) - 180
        lref = lon_0[ np.shape(lon_0)[0]/2 ]
        lon_0 = np.mod(lon_0 - (lref - 180), 360) + (lref - 180)
        #lon_0 = np.mod(lon_0 - np.min(lon_0), 360) - lref # - 180
        lon_0 = np.rad2deg(np.unwrap(np.deg2rad(lon_0)))
        lat_0 = lat_all[ind_ipass[i]:ind_ipass[i + 1]]
        var_0 = var_all[ind_ipass[i]:ind_ipass[i + 1]]
        time_0 = time_all[ind_ipass[i]:ind_ipass[i + 1]]
        ipass_0 = ipass_all[ind_ipass[i]:ind_ipass[i + 1]]
        cycle_0 = cycle_all[ind_ipass[i]:ind_ipass[i + 1]]
        dtime = time_0[1:] - time_0[:-1]
        dlon = abs(lon_0[1:] - lon_0[:-1])
        if (np.array(dlon) > 180).any():
            lon_0 = np.mod(lon_0, 360)
        if len(time_0) > 3:
            delta = stats.mode(dtime)[0][0]
        else:
            continue
        ndelta = np.round(dtime / delta).astype('int')
        ind_dtime = np.where(ndelta >= 2)[0]
        if ind_dtime.size != 0:
            time = time_0[:ind_dtime[0] + 1]
            var = var_0[:ind_dtime[0] + 1]
        else:
            time = time_0
            var = var_0
        for i in range(len(ind_dtime)):
            # time_fill = np.linspace(time_0[ind_dtime[i]],
            #                         time_0[ind_dtime[i] + 1],
            #                         num=(time_0[ind_dtime[i] + 1]
            #                         - time_0[ind_dtime[i]]) / delta)
            time_fill = np.linspace(time_0[ind_dtime[i]],
                                    time_0[ind_dtime[i] + 1],
                                    num=ndelta[ind_dtime[i]], endpoint=False)
            var_fill = np.zeros(np.shape(time_fill[1:])) * np.nan
            # time = time.append(time_fill[1:])
            time = np.hstack([time, time_fill[1:]])
            # var = var.append(var_fill)
            var = np.hstack([var, var_fill])
            if i != (len(ind_dtime) - 1):
                time = np.hstack([time, time_0[ind_dtime[i]
                                  + 1:ind_dtime[i + 1] + 1]])
                var = np.hstack([var, var_0[ind_dtime[i]
                                 + 1:ind_dtime[i + 1] + 1]])
            else:
                time = np.hstack([time, time_0[ind_dtime[i] + 1:]])
                var = np.hstack([var, var_0[ind_dtime[i] + 1:]])
        # time = np.concatenate(time)
        # var = np.concatenate(var)
        func = interpolate.interp1d(time_0, lon_0, kind='quadratic')
        lon = func(time)
        func = interpolate.interp1d(time_0, lat_0, kind='quadratic')
        lat = func(time)
        ssha = var
        mask_gap = np.isnan(ssha)
        lon = lon - np.floor((np.min(lon) + 180.) / 360.) * 360.
        time = np.float64(time)
        start_time = num2date(time[0], time_units)
        end_time = num2date(time[-1], time_units)
        time = num2date(time, units=time_units)
        time_num = time
        for t in range(np.shape(time)[0]):
            time_num[t] = date2num(time[t], units='microseconds since 1970-01-01 00:00:00.000000Z')
            time_num[t] *= 10**(-6)
        #time = (num2date(time, time_units)
        #        - calendar.timegm(cnes_day.timetuple()))
        # for t in range(np.shape(time)[0]):
        #    time[t] = date2num(num2date(time[t], units) - cnes_day, unit=)
        #time = time.total_seconds()
        #time = (time + calendar.timegm(julian_day.timetuple())
        #        - calendar.timegm(cnes_day.timetuple()))
        # NOTE : lon/lat must be continuous even if crossing dateline
        # (ie. no [-180,180] clipping)
        # Make GCPs (mimic a swath of arbitrary width in lon/lat, here ~5km)
        # gcps = tools_for_gcp.make_gcps_v1(lon, lat, dist_gcp=dist_gcp)
        gcps = tools_for_gcp.make_gcps_v2(lon, lat, dist_gcp=dist_gcp)
        # Write geotiff
        # NOTE : product_name to be changed, set here for test
        metadata = {}
        (dtime, time_range) = stfmt.format_time_and_range(start_time,
                                                          end_time, units='s')
        metadata['product_name'] = L3_MAPS[L3id]['productname']
        dname = (os.path.splitext(os.path.basename(infile))[0]
                 + '_c' + str(int(cycle_0[0])).zfill(4)
                 + '_p' + str(int(ipass_0[0])).zfill(3))
        metadata['name'] = dname
        metadata['datetime'] = dtime
        metadata['time_range'] = time_range
        metadata['begin_time'] =  start_time.strftime(TIMEFMT)
        metadata['end_time'] = end_time.strftime(TIMEFMT)
        metadata['source_URI'] = infile
        metadata['source_provider'] = 'AVISO'
        metadata['processing_center'] = ''
        metadata['conversion_version'] = '0.0.0'
        metadata['conversion_datatime'] = stfmt.format_time(datetime.utcnow())
        metadata['type'] = 'along_track'
        metadata['cycle'] = int(cycle_0[0])
        metadata['pass'] = int(ipass_0[0])
        metadata['spatial_resolution'] = np.float32(7000)
        geolocation = {}
        geolocation['projection'] = stfmt.format_gdalprojection()
        geolocation['gcps'] = stfmt.format_gdalgcps(*gcps)
        band = []
        mask = np.ma.getmaskarray(ssha)
        if write_netcdf:
            vmin = np.nanmin(ssha)
            vmin_pal = vmin
            vmax = np.nanmax(ssha)
            vmax_pal = vmax
            #print('bla')
        scale = (vmax - vmin) / 254.
        offset = vmin
        mask = np.ma.getmaskarray(ssha)
        array = np.clip(np.round((ssha - offset) / scale),
                        0, 254).astype('uint8')
        array[mask] = 255
        array[mask_gap] = 255
        if write_netcdf is False:
            array = array[:, np.newaxis]
        colortable = stfmt.format_colortable('matplotlib_jet',
                                             vmax=vmax, vmax_pal=vmax_pal,
                                             vmin=vmin, vmin_pal=vmin_pal)
        band.append({'array': array, 'scale': scale, 'offset': offset,
                     'description': L3_MAPS[L3id]['parameter'],
                     'name': L3_MAPS[L3id]['hname'],
                     'unittype': 'm', 'nodatavalue': 255,
                     'parameter_range': [vmin, vmax],
                     'colortable': colortable})
        if write_netcdf is False:
            tifffile = stfmt.format_tifffilename(outdir, metadata,
                                                 create_dir=True)
            stfmt.write_geotiff(tifffile, metadata, geolocation, band)
        else:
           geolocation['geotransform'] = [lon, lat]
           geolocation['time'] = time_num[:]
           netcdffile = stfmt.format_ncfilename(outdir, metadata,
                                                create_dir=True)
           write_netcdf_1d(netcdffile, metadata, geolocation, band,
                          model='along_track', dgcps=1)
        # Check GCPs
        # tools_for_gcp.check_gcps(tifffile, lon, lat, gcps[0], gcps[1])
        # if ipass_0[0] == 672 or ipass_0[0] == 676:
        #     print tifffile
        #     tools_for_gcp.check_gcps(tifffile, lon, lat, gcps[0], gcps[1])

