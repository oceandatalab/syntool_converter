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

from cerbere.mapper.ghrsstncfile import GHRSSTNCFile
import numpy as np
from netCDF4 import num2date
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os


def ghrsst_oi_mw_sst(infile, outdir,
                vmin=271.05, vmax=309.15, vmin_pal=273., vmax_pal=305.):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    ghrsst = GHRSSTNCFile(infile)
    sst = ghrsst.read_values('analysed_sst')[::-1, :]
    lon = ghrsst.read_values('lon')[:]
    indsorted = np.argsort(lon)
    sst = sst[:,indsorted]
    mask = ghrsst.read_values('mask')[::-1, :]
    mask = mask[:,indsorted]
    #sea_ice_fraction = ghrsst.read_values('sea_ice_fraction')[::-1, :]
    # lon = ghrsst.read_values('lon')
    # dlon = lon[1] - lon[0]
    # lon0 = lon[0] - dlon / 2
    # lat = ghrsst.read_values('lat')[::-1]
    # dlat = lat[1] - lat[0]
    # lat0 = lat[0] - dlat / 2
    lon0 = ghrsst.read_global_attribute('westernmost_longitude')
    dlonstring = ghrsst.read_global_attribute('geospatial_lon_resolution')
    dlon = float(dlonstring.strip(" deg"))
    lat0 = ghrsst.read_global_attribute('northernmost_latitude')
    dlatstring = ghrsst.read_global_attribute('geospatial_lat_resolution')
    dlat = -float(dlatstring.strip(" deg"))
    dtime = ghrsst.read_values('time')[0]
    dtime_units = ghrsst.read_field('time').units
    dtime = num2date(dtime, dtime_units)
    ghrsst_id = ghrsst.read_global_attribute('id')
    if 'glob' in ghrsst_id.lower():
        product_name = 'REMSS_MWOI_SST'
    else:
        raise Exception('Unknown REMSS ID : {}'.format(ghrsst_id))

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = product_name
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = ['-12h', '+12h']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'Ifremer'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea surface temperature'
    metadata['type'] = 'remote sensing'
    metadata['longitude_resolution'] = abs(dlon)
    metadata['latitude_resolution'] = abs(dlat)
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon0-dlon/2, dlon, 0,
                                   lat0-dlat/2, 0, dlat]
    band = []
    #indndv = np.where((sst.mask == True) | (sea_ice_fraction > 0))
    indndv = np.where((sst.mask == True) | (mask != 65))
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(sst, vmin, vmax, out=sst)
    array = np.round((sst - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('cerbere_medspiration',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea surface temperature', 'unittype':'K',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
