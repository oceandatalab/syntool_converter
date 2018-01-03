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

from cerbere.mapper.ncfile import NCFile
import numpy as np
from netCDF4 import num2date
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os

def bathymetry_gebco(infile, outdir,
                vmin=-6000, vmax=0, vmin_pal=-6000., vmax_pal=0.):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    ncfile = NCFile(infile)
    bat = ncfile.read_values('elevation')[:, :]
    bat=bat.astype('float32')
    bat[(bat<0) & (bat>=-25)]=-25
    bat[(bat<-25) & (bat>=-50)]=-50
    bat[(bat<-50) & (bat>=-100)]=-100
    bat[(bat<-100) & (bat>=-500)]=-500
    bat[(bat<-500) & (bat>=-1000)]=-1000
    bat[(bat<-1000) & (bat>=-2000)]=-2000
    bat[(bat<-2000) & (bat>=-3000)]=-3000
    bat[(bat<-3000) & (bat>=-4000)]=-4000
    bat[(bat<-4000) & (bat>=-5000)]=-5000
    bat[(bat<-5000) & (bat>=-6000)]=-6000
    bat[(bat<-6000) & (bat>=-10000)]=-10000
    mask = [bat>=0]
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(bat, vmin, vmax, out=bat)
    bat-=offset
    bat/=scale
    bat = np.round(bat).astype('uint8')
    lon = ncfile.read_values('lon')[:2:1]
    lat = ncfile.read_values('lat')[:2:1]
    lon0 = lon[0]
    dlon = lon[-1]-lon[0]
    lat0 = lat[0]
    dlat = lat[-1]-lat[0]
   
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = 'GEBCO bathymetry'
    metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['datetime'] = stfmt.format_time(datetime(2012,1,1))
    metadata['time_range'] = ['-3660d', '+3660d']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'GEBCO'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'bathymetry '
    metadata['type'] = 'remote sensing'
    metadata['longitude_resolution'] = abs(dlon)
    metadata['latitude_resolution'] = abs(dlat)
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon0-dlon/2., dlon, 0,
                                   lat0-dlat/2., 0, dlat]
    band = []
    bat[mask] = 255
    colortable = stfmt.format_colortable('ibcso',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)

    band.append({'array':bat, 'scale':scale, 'offset':offset,
                 'description':'bathymetry', 'unittype':'m',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
