# -*- coding: utf-8 -*-

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

import netCDF4
import numpy
import os
import syntool_converter.utils.syntoolformat as stfmt
import datetime
import logging

logger = logging.getLogger(__name__)


def smos_l3_bec_sss(infile, outdir,
                vmin=31.825, vmax=38.175, vmin_pal=32, vmax_pal=38):
    """
    """
    # Read/Process data
    logger.info('Read/Process data')
    smos = netCDF4.Dataset(infile, 'r')
    time_start = netCDF4.num2date(smos['time'][0], smos['time'].units)
    time_stop = time_start + datetime.timedelta(days=1)
    lat = smos['lat'][::-1]
    lon = smos['lon'][:]
    sss = smos['oa_sss'][0, ::-1, :]

    # Construct metadata/geolocation/band(s)
    logger.info('Construct metadata/geolocation/band(s)')
    dtime, time_range = stfmt.format_time_and_range(time_start, time_stop,
                                                    units='h')
    lat0, dlat = lat[0], lat[1] - lat[0]
    lon0, dlon = lon[0], lon[1] - lon[0]
    now = datetime.datetime.utcnow()
    metadata = {}
    metadata['product_name'] = 'SMOS_L3_BEC_SSS'
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    src = 'SMOS Barcelona Expert Centre, ICM-CSIC / UPC, Barcelona, Spain'
    metadata['source_provider'] = src
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(now)
    metadata['parameter'] = 'sea surface salinity'
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon0-dlon/2., dlon, 0,
                                   lat0-dlat/2., 0, dlat]
    band = []
    offset, scale = vmin, (vmax-vmin)/254.
    numpy.clip(sss.data, vmin, vmax, out=sss.data)
    array = numpy.round((sss.data - offset) / scale).astype('uint8')
    array[sss.mask] = 255
    colortable = stfmt.format_colortable('matplotlib_jet',
                                         vmin=vmin, vmax=vmax,
                                         vmin_pal=vmin_pal, vmax_pal=vmax_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea surface salinity', 'unittype':'PSS',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
    logger.info('Write geotiff')
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
