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

from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import syntool_converter.utils.syntoolformat as stfmt
import os


def smosstorm_smap_wind(infile, outdir,
                        vmin=0., vmax=50.8, vmin_pal=0., vmax_pal=50.):
    """
    """
    # Read data
    dataset = Dataset(infile)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    minute = dataset.variables['minute'][:]
    ## for wind we disable auto masking because of wrong valid_min/valid_max
    ## (valid_min/valid_max should be set for values before offset/scale)
    dataset.variables['wind'].set_auto_mask(False)
    wind = dataset.variables['wind'][:]
    wind = np.ma.masked_less(wind, 0)
    day = datetime(dataset.year_of_observation,
                   dataset.month_of_observation,
                   dataset.getncattr('day_of_month_of observation'))
    dataset.close()
    name_pattern = '{}_{{}}_{{}}'.format(os.path.splitext(os.path.basename(infile))[0])
    if np.ma.is_masked(lat):
        raise Exception('Some lat are masked.')
    dlat = lat[1] - lat[0]
    if not np.all(np.isclose(lat[1:] - lat[:-1], dlat)):
        raise Exception('Unexpected not unique dlat.')
    if np.ma.is_masked(lon):
        raise Exception('Some lon are masked.')
    dlon = lon[1] - lon[0]
    if not np.all(np.isclose(lon[1:] - lon[:-1], dlon)):
        raise Exception('Unexpected not unique dlon.')
    if np.abs(lon[-1] + dlon - lon[0]) != 360:
        raise Exception('lon does not cover exactly 360 degrees.')
    if np.all(np.ma.getmaskarray(minute) | np.ma.getmaskarray(wind)):
        raise Exception('minute + wind is fully masked.')
    if minute.min() < 0 or minute.max() > 1440:
        raise Exception('Unexpected minute outside [0, 1440].')

    # Flip grid
    if dlat > 0:
        lat = lat[::-1]
        dlat = -dlat
        minute = minute[::-1, :, :]
        wind = wind[::-1, :, :]
    if dlon < 0:
        lon = lon[::-1]
        dlon = -dlon
        minute = minute[:, ::-1, :]
        wind = wind[:, ::-1, :]

    # Loop on node (ascending / descending) and pass
    nlon = lon.size
    pass_dminute = 49
    for inode in [0, 1]:
        valid_minute = minute[:, :, inode].compressed()
        if valid_minute.size == 0:
            continue
        hist, _ = np.histogram(valid_minute, bins=1441, range=(0, 1441))
        indminute = np.where(hist != 0)[0]
        dminute = indminute[1:] - indminute[:-1]
        indsplit = np.where(dminute > pass_dminute)[0]
        minute0 = np.concatenate((indminute[[0]], indminute[indsplit + 1]))
        minute1 = np.concatenate((indminute[indsplit], indminute[[-1]]))
        for _minute0, _minute1 in zip(minute0, minute1):
            indpass = np.where((minute[:, :, inode] >= _minute0) & \
                               (minute[:, :, inode] <= _minute1) & \
                               (np.ma.getmaskarray(wind[:, :, inode]) == False))
            if indpass[0].size == 0:
                continue
            lat_slice = slice(indpass[0].min(), indpass[0].max() + 1)
            indlon = np.unique(indpass[1])
            if indlon.size == 1:
                dindlon = np.array([0])
            else:
                dindlon = indlon[1:] - indlon[:-1]
            if (indlon[-1] - indlon[0]) <= (nlon - dindlon.max()):
                lon_slice = slice(indlon[0], indlon[-1] + 1)
                lat0 = lat[lat_slice.start]
                lon0 = lon[lon_slice.start]
                _wind = wind[lat_slice, lon_slice, inode].copy()
                _minute = minute[lat_slice, lon_slice, inode].copy()
            else:
                indsplit = dindlon.argmax()
                lon_slice0 = slice(indlon[indsplit + 1], nlon)
                lon_slice1 = slice(0, indlon[indsplit] + 1)
                lat0 = lat[lat_slice.start]
                lon0 = lon[lon_slice0.start]
                _wind = np.ma.concatenate((wind[lat_slice, lon_slice0, inode],
                                           wind[lat_slice, lon_slice1, inode]),
                                          axis=1)
                _minute = np.ma.concatenate((minute[lat_slice, lon_slice0, inode],
                                             minute[lat_slice, lon_slice1, inode]),
                                            axis=1)
            lon0 = np.mod(lon0 + 180., 360.) - 180.
            _wind[_minute < _minute0] = np.ma.masked
            _wind[_minute > _minute1] = np.ma.masked
            start_time = day + timedelta(seconds=_minute0 * 60)
            stop_time = day + timedelta(seconds=_minute1 * 60)
            name = name_pattern.format(start_time.strftime('%H%M'),
                                       stop_time.strftime('%H%M'))

            # Construct metadata/geolocation/band(s)
            print 'Construct metadata/geolocation/band(s)'
            metadata = {}
            (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                              units='ms')
            metadata['product_name'] = 'SMOSSTORM_SMAP_wind'
            metadata['name'] = name
            metadata['datetime'] = dtime
            metadata['time_range'] = time_range
            metadata['source_URI'] = infile
            metadata['source_provider'] = 'NASA'
            metadata['processing_center'] = 'Remote Sensing Systems'
            metadata['conversion_software'] = 'Syntool'
            metadata['conversion_version'] = '0.0.0'
            metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
            metadata['parameter'] = 'wind speed'
            geolocation = {}
            geolocation['projection'] = stfmt.format_gdalprojection()
            geolocation['geotransform'] = [lon0 - dlon / 2., dlon, 0,
                                           lat0 - dlat / 2., 0, dlat]
            band = []
            mask = np.ma.getmaskarray(_wind)
            _wind = np.ma.getdata(_wind)
            offset, scale = vmin, (vmax - vmin) / 254.
            np.clip(_wind, vmin, vmax, out=_wind)
            array = np.round((_wind - offset) / scale).astype('uint8')
            array[mask] = 255
            colortable = stfmt.format_colortable('matplotlib_jet',
                                                 vmax=vmax, vmax_pal=vmax_pal,
                                                 vmin=vmin, vmin_pal=vmin_pal)
            band.append({'array':array, 'scale':scale, 'offset':offset,
                         'description':'wind speed', 'unittype':'m/s',
                         'nodatavalue':255, 'parameter_range':[vmin, vmax],
                         'colortable':colortable})

            # Write geotiff
            print 'Write geotiff'
            tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
            stfmt.write_geotiff(tifffile, metadata, geolocation, band)
