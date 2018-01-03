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
import numpy.ma as ma
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os


def smosstorm_ascat_wind(infile, outdir,
                         vmin=0., vmax=50.8, vmin_pal=0., vmax_pal=50.):
    """
    """
    # Read/Process data
    dataset = Dataset(infile)
    lon = dataset.variables['lon'][:]
    valid = np.where(ma.getmaskarray(lon) == False)
    slices = [slice(valid[0].min(), valid[0].max() + 1),
              slice(valid[1].min(), valid[1].max() + 1)]
    lon = lon[slices]
    if lon.shape[1] % 2 != 0:
        raise Exception('Number of cells should be even.')
    swath_ncell = lon.shape[1] / 2
    lat = dataset.variables['lat'][slices]
    wind_speed = dataset.variables['wind_speed'][slices]
    wind_dir = dataset.variables['wind_dir'][slices]
    time = dataset.variables['time'][slices]
    time_units = dataset.variables['time'].units
    start_time = num2date(time.min(), time_units)
    stop_time = num2date(time.max(), time_units)
    datagroup = '_'.join(os.path.splitext(os.path.basename(infile))[0].split('_')[0:3])
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SMOSSTORM_ASCAT_wind'
    metadata['datagroup'] = datagroup
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'EUMETSAT'
    metadata['processing_center'] = 'KNMI'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = ['wind speed', 'wind direction']
    for i in range(2):
        if i == 0:
            metadata['name'] = datagroup + '_left'
        else:
            metadata['name'] = datagroup + '_right'
        swath_slice = slice(i * swath_ncell, (i + 1) * swath_ncell)
        swath_lon = lon[:, swath_slice]
        swath_lat = lat[:, swath_slice]
        swath_wind_speed = wind_speed[:, swath_slice]
        swath_wind_dir = wind_dir[:, swath_slice]
        dgcp = 16. # spacing is about 12.5km
        ngcps = np.ceil(np.array(swath_lon.shape) / dgcp) + 1.
        pix = np.linspace(0, swath_lon.shape[1] - 1, num=ngcps[1]).round().astype('int32')
        lin = np.linspace(0, swath_lon.shape[0] - 1, num=ngcps[0]).round().astype('int32')
        pix2d, lin2d = np.meshgrid(pix, lin)
        gcplon = swath_lon[lin2d, pix2d]
        gcplon[np.where(gcplon >= 180)] -= 360
        if gcplon.max() - gcplon.min() > 180:
            gcplon[np.where(gcplon < 0)] += 360.
        gcplat = swath_lat[lin2d, pix2d]
        gcppix = pix2d + 0.5
        gcplin = lin2d + 0.5
        gcphei = np.zeros(ngcps)
        geolocation = {}
        geolocation['projection'] = stfmt.format_gdalprojection()
        geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                    gcppix, gcplin)
        band = []
        indndv = np.where(ma.getmaskarray(swath_wind_speed))
        offset, scale = vmin, (vmax-vmin)/254.
        clipped = np.clip(ma.getdata(swath_wind_speed), vmin, vmax)
        array = np.round((clipped - offset) / scale).astype('uint8')
        array[indndv] = 255
        colortable = stfmt.format_colortable('matplotlib_jet',
                                             vmax=vmax, vmax_pal=vmax_pal,
                                             vmin=vmin, vmin_pal=vmin_pal)
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'wind speed', 'unittype':'m/s',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax],
                     'colortable':colortable})
        indndv = np.where(ma.getmaskarray(swath_wind_dir))
        clipped = np.clip(ma.getdata(swath_wind_dir), 0, 360)
        array = np.round(clipped / 360. * 254.).astype('uint8')
        array[indndv] = 255
        band.append({'array':array, 'scale':360./254., 'offset':0.,
                     'description':'wind direction', 'unittype':'deg',
                     'nodatavalue':255, 'parameter_range':[0, 360.]})
        # Write geotiff
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
