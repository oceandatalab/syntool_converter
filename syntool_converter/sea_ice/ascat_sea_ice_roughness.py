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
from datetime import datetime
import numpy as np
import syntool_converter.utils.syntoolformat as stfmt
import os
import osr


# def test_xy(dset, proj, tf):
#     """
#     """
#     import pyproj
#     import matplotlib.pyplot as plt
#     import pdb
#     lon = dset.variables['longitude'][0, :, :]
#     lat = dset.variables['latitude'][0, :, :]
#     proj = pyproj.Proj(proj)
#     x, y = proj(lon, lat)
#     xproj = np.arange(x.shape[1]) * tf[1] + tf[0] + tf[1] / 2.
#     yproj = np.arange(x.shape[0]) * tf[5] + tf[3] + tf[5] / 2.
#     plt.figure()
#     plt.imshow(np.sqrt((y - yproj[:, np.newaxis]) ** 2. + (x - xproj[np.newaxis, :]) ** 2.))
#     plt.colorbar()
#     #plt.show()
#     #pdb.set_trace()
#     x0 = x[:,0].mean() - 12500. / 2.
#     x1 = x[:,-1].mean() + 12500. / 2.
#     xd = x1 - x0 - 12500 * x.shape[1]
#     y0 = y[0,:].mean() + 12500. / 2.
#     y1 = y[-1,:].mean() - 12500. / 2.
#     yd = y0 - y1 - 12500 * y.shape[0]
#     tf[0] = x0 + xd / 2.
#     tf[3] = y0 - yd / 2.
#     xproj = np.arange(x.shape[1]) * tf[1] + tf[0] + tf[1] / 2.
#     yproj = np.arange(x.shape[0]) * tf[5] + tf[3] + tf[5] / 2.
#     plt.figure()
#     plt.imshow(np.sqrt((y - yproj[:, np.newaxis]) ** 2. + (x - xproj[np.newaxis, :]) ** 2.))
#     plt.colorbar()
#     plt.show()
#     pdb.set_trace()


def ascat_sea_ice_roughness(infile, outdir,
                            vmin=0., vmax=0.122409712058,
                            vmin_pal=0., vmax_pal=0.122409712058):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    dset = Dataset(infile)
    dset.variables['sigma_40_mask'].set_auto_maskandscale(False)
    roughness = dset.variables['sigma_40_mask'][0, :, :].astype('uint8')
    indnodata = np.where((roughness == 0) | (roughness == 255))
    roughness = roughness * dset.variables['sigma_40_mask'].scale_factor + \
                dset.variables['sigma_40_mask'].add_offset
    dtime = num2date(dset.variables['time'][0], dset.variables['time'].units)
    pole = dset.pole
    name = os.path.splitext(os.path.basename(infile))[0]
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = 'ASCAT_sea_ice_roughness'
    metadata['name'] = name
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = ['-12h', '+12h']
    metadata['source_URI'] = infile
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea ice roughness'
    geolocation = {}
    srs = osr.SpatialReference()
    if pole == 'north':
        srs.ImportFromEPSG(3411) # use 3413 if problems with ellipsoid
        geolocation['projection'] = srs.ExportToWkt()
        geolocation['geotransform'] = [-3850000, 12500, 0, 5850000, 0, -12500]
    elif pole == 'south':
        srs.ImportFromEPSG(3412) # use 3976 if problems with ellipsoid
        geolocation['projection'] = srs.ExportToWkt()
        geolocation['geotransform'] = [-3950000, 12500, 0, 4350000, 0, -12500]
    else:
        raise Exception('Which pole ?')
    #test_xy(dset, srs.ExportToProj4(), geolocation['geotransform'])
    band = []
    np.clip(roughness, vmin, vmax, out=roughness)
    offset, scale = vmin, (vmax - vmin) / 254.
    array = np.round((roughness - offset) / scale).astype('uint8')
    array[indnodata] = 255
    colortable = stfmt.format_colortable('pesket', vmin=vmin, vmax=vmax,
                                         vmin_pal=vmin_pal, vmax_pal=vmax_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea ice roughness', 'unittype':'',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
