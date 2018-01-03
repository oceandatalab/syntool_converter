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

from pyhdf.SD import SD, SDC
from datetime import datetime, timedelta
import numpy as np
import syntool_converter.utils.syntoolformat as stfmt
import os
import osr
#import pdb


# File patterns :
# asi-s6250-20020808-v5.hdf
# asi-AMSR2-s6250-20140922-v5.hdf
# asi-n6250-20020601-v5.hdf
# asi-AMSR2-n6250-20140922-v5.hdf


def amsr_sea_ice_conc(infile, outdir,
                      vmin=0., vmax=100., vmin_pal=0., vmax_pal=100.,
                      write_netcdf=False):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    hsd = SD(infile, SDC.READ)
    if 'ASI Ice Concentration (Version: 5.3)' in hsd.datasets():
        sds = hsd.select('ASI Ice Concentration (Version: 5.3)')
    elif 'ASI Ice Concentration':
        sds = hsd.select('ASI Ice Concentration')
    else:
        raise Exception('No supported dataset available')
    attrs = sds.attributes()
    long_name = attrs['long_name']
    long_name_split = long_name.split(', ')

    date_str = ''
    sensor_name = ''
    if 5 == len(long_name_split):
        date_str = long_name_split[1]
        sensor_name = long_name_split[3]
    elif 6 == len(long_name_split):
        date_str = long_name_split[2]
        sensor_name = long_name_split[4]
    else:
        raise Exception('long_name does not have the expected format.')

    dtime = datetime.strptime(date_str, '%Y%m%d') + timedelta(0.5)
    name = os.path.splitext(os.path.basename(infile))[0]
    if sensor_name == 'AMSR-E':
        name = name.replace('asi', 'asi-AMSRE') # make uniform with AMSR2
        sensor_platform = 'Aqua'
    elif sensor_name == 'AMSR2':
        sensor_platform = 'GCOM-W1'
    else:
        raise Exception('Unexpected sensor.')
    pole = name.split('-')[2][0]
    conc = sds.get()[::-1, :]
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = 'AMSR_sea_ice_concentration'
    metadata['name'] = name
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = ['-12h', '+12h']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'Universität Bremen'
    metadata['processing_center'] = 'Universität Bremen'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea ice concentration'
    metadata['type'] = 'remote sensing'
    metadata['sensor_type'] = 'radiometer'
    metadata['sensor_name'] = sensor_name
    metadata['sensor_platform'] = sensor_platform
    geolocation = {}
    srs = osr.SpatialReference()
    if pole == 'n':
        srs.ImportFromEPSG(3411) # use 3413 if problems with ellipsoid
        geolocation['projection'] = srs.ExportToWkt()
        geolocation['geotransform'] = [-3850000, 6250, 0, 5850000, 0, -6250]
    elif pole == 's':
        srs.ImportFromEPSG(3412) # use 3976 if problems with ellipsoid
        geolocation['projection'] = srs.ExportToWkt()
        geolocation['geotransform'] = [-3950000, 6250, 0, 4350000, 0, -6250]
    else:
        raise Exception('Which pole ?')
    band = []
    conc[np.where(~np.isfinite(conc))] = -1 # land
    conc[np.where(conc == 0)] = -1 # no ice at all
    ndv = np.where(conc < vmin)
    np.clip(conc, vmin, vmax, out=conc)
    offset, scale = vmin, (vmax-vmin)/254.
    array = np.round((conc - offset) / scale).astype('uint8')
    array[ndv] = 255
    colortable = stfmt.format_colortable('sea_ice_conc', vmin=vmin, vmax=vmax,
                                         vmin_pal=vmin_pal, vmax_pal=vmax_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea ice concentration', 'unittype':'%',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    # Write
    if write_netcdf == False:
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    elif write_netcdf == True:
        print 'Write netcdf'
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        band[0]['name'] = 'sea_ice_concentration'
        band[0]['long_name'] = 'ASI Ice Concentration'
        band[0]['unittype'] = 'percent'
        metadata['spatial_resolution'] = 6250.
        dgcps = np.round(100. / np.array([6.250, 6.250])).astype('int')
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'grid_proj',
                           dgcps=dgcps)
