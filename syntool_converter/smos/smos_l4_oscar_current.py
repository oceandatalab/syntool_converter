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

from cerbere.mapper.smosncfile import SMOSNCFile
from netCDF4 import num2date
import numpy as np
import os
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime, timedelta


def smos_l4_oscar_current(infile, outdir,
                          vmin=0., vmax=5.08, vmin_pal=0., vmax_pal=2.):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    smos = SMOSNCFile(infile)
    time_start_units = smos.read_field('date_start').units
    time_start = num2date(smos.read_values('date_start')[0], time_start_units)
    time_stop_units = smos.read_field('date_stop').units
    time_stop = num2date(smos.read_values('date_stop')[0], time_stop_units)
    if (time_stop - time_start) == timedelta(days=6):
        time_stop += timedelta(days=1)
    lat = smos.read_values('lat')[::-1]
    lon = smos.read_values('lon')
    ucur = smos.read_values('Zonal_component_surface_currents')[::-1, :]
    vcur = smos.read_values('Meridional_component_surface_currents')[::-1, :]

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    dtime, time_range = stfmt.format_time_and_range(time_start, time_stop,
                                                    units='h')
    lat0, dlat = lat[0], lat[1] - lat[0]
    lon0, dlon = lon[0], lon[1] - lon[0]
    metadata = {}
    metadata['product_name'] = 'SMOS_L4_OSCAR_current'
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'Ifremer/CNES'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = ['current velocity', 'current direction']
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
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
