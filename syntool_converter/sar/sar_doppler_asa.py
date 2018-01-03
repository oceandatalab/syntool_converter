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
from datetime import datetime, timedelta
import syntool_converter.utils.syntoolformat as stfmt
from matplotlib.colors import LinearSegmentedColormap, Normalize


# def test_gdalgcps():
#     infile = '/local/home/fab/data/sar/ASA/agulhas/ASA_WSM_1PNPDE20110824_211403_000002143106_00014_49597_2093/SAR_doppler.nc'
#     sardop = NCFile(infile)
#     lon = sardop.read_values('longitude')[::-1, :]
#     lat = sardop.read_values('latitude')[::-1, :]
#     shp = lon.shape
#     import gdal
#     #infile = '/local/home/data/tmp/ASA_WSM_1PNPDE20110824_211403_000002143106_00014_49597_2093/sar_doppler.tiff'
#     infile = '/local/home/data/syntool_inputs/sar_doppler/2011/08/24/ASA_WSM_1PNPDE20110824_211403_000002143106_00014_49597_2093/sar_doppler.tiff'
#     src_ds, dst_ds = gdal.Open(infile), None
#     #options = ['']
#     #options = ['MAX_GCP_ORDER=3']
#     options = ['MAX_GCP_ORDER=-1']
#     transformer = gdal.Transformer(src_ds, dst_ds, options)
#     xy_grid = np.mgrid[slice(0, shp[1]), slice(0, shp[0])].astype('int32')
#     xy_dims = xy_grid.shape[1:3]
#     xy_grid = xy_grid.reshape((2, -1)).transpose()
#     lonlat = transformer.TransformPoints(0, xy_grid)
#     lonlat = np.array(lonlat[0], dtype='float32')
#     longdal = lonlat[:, 0].reshape(xy_dims).transpose()
#     latgdal = lonlat[:, 1].reshape(xy_dims).transpose()
#     import matplotlib.pyplot as plt
#     plt.plot(lon.flatten(), lat.flatten(), '+b')
#     plt.plot(longdal.flatten(), latgdal.flatten(), '+r')
#     plt.show()
#     plt.close()
#     import pdb
#     pdb.set_trace()


def doppler_colormap():
    """
    """
    xval = [-2.5, -2.0, -1.6, -0.8, -0.3, 0, 0.3, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5]
    red = [255, 220, 128, 12, 50, 255, 90, 255, 255, 255, 255, 255, 160]
    green = [0, 0, 30, 120, 200, 255, 255, 255, 255, 170, 85, 0, 0]
    blue = [128, 220, 200, 238, 235, 255, 30, 0, 0, 0, 0, 0, 0]
    segmentdata = {'red':[], 'green':[], 'blue':[]}
    for x, r, g, b in zip(xval, red, green, blue):
        xnorm = (x+2.5)/5
        rnorm = r/255.
        gnorm = g/255.
        bnorm = b/255.
        segmentdata['red'].append((xnorm, rnorm, rnorm))
        segmentdata['green'].append((xnorm, gnorm, gnorm))
        segmentdata['blue'].append((xnorm, bnorm, bnorm))
    return LinearSegmentedColormap('doppler', segmentdata)


def sar_doppler(infile, outdir):
    """
    """
    # tmp
    #infile = '/local/home/fab/data/sar/ASA/agulhas/ASA_WSM_1PNPDE20110518_210602_000002143102_00330_48189_1274/SAR_doppler.nc'
    # infile = '/local/home/fab/data/sar/ASA/agulhas/ASA_WSM_1PNPDE20110824_211403_000002143106_00014_49597_2093/SAR_doppler.nc'
    # outdir = '/local/home/data/syntool_inputs'
    # /tmp
    # Read/Process data
    print 'Read/Process data'
    sardop = NCFile(infile)
    product_ref = sardop.read_global_attribute('SOURCE_PRODUCT_REF')
    start_time = sardop.read_global_attribute('SOURCE_START_DATE')
    start_time = datetime.strptime(start_time, '%Y%m%d%H%M%S.%f')
    duration = sardop.read_global_attribute('SOURCE_ACQ_DURATION')
    stop_time = start_time + timedelta(seconds=duration)
    polarisation = sardop.read_global_attribute('SOURCE_POLARIZATION')
    lon = sardop.read_values('longitude')[::-1, :]
    lat = sardop.read_values('latitude')[::-1, :]
    #dopano = sardop.read_values('dopanomaly')[::-1, :]
    radvel = sardop.read_values('radial_vel')[::-1, :]
    validity = sardop.read_values('validity')[::-1, :]
    track_angle = sardop.read_global_attribute('SOURCE_TRACK_ANGLE')
    if track_angle < 0:
        radvel *= -1
    shp = lon.shape
    nlines = np.ceil(shp[0]/4.)+1
    lines = np.round(np.linspace(0, shp[0]-1, num=nlines)).astype('int32')
    npixels = np.ceil(shp[1]/4)+1
    pixels = np.round(np.linspace(0, shp[1]-1, num=npixels)).astype('int32')
    gcplin = np.tile(lines.reshape(nlines, 1), (1, npixels))
    gcppix = np.tile(pixels.reshape(1, npixels), (nlines, 1))
    gcplon = lon[gcplin, gcppix]
    gcplat = lat[gcplin, gcppix]
    gcphei = np.zeros((nlines, npixels))
    gcppix = gcppix + 0.5
    gcplin = gcplin + 0.5
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                        units='ms')
    metadata['product_name'] = 'SAR_doppler'
    metadata['name'] = product_ref
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'ESA'
    metadata['processing_center'] = 'CLS'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'radial horizontal velocities'
    metadata['type'] = 'remote sensing'
    metadata['sensor_type'] = 'SAR'
    metadata['sensor_name'] = 'ASAR'
    metadata['sensor_platform'] = 'ENVISAT'
    metadata['sensor_mode'] = 'WSM'
    #metadata['sensor_swath'] = sensor_swath
    metadata['sensor_polarisation'] = polarisation
    #metadata['sensor_pass'] = sensor_pass
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection(geogcs='WGS84')
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                gcppix, gcplin)
    # band = []
    # scale = (vmax-vmin)/254.
    # offset = vmin
    # indzero = np.where(validity == 0)
    # array = np.clip(np.round((radvel-offset)/scale), 0, 254).astype('uint8')
    # array[indzero] = 255
    # band.append({'array':array, 'scale':scale, 'offset':offset,
    #              'description':'radial horizontal velocities', 'unittype':'m/s',
    #              'nodatavalue':255, 'parameter_range':[vmin, vmax]})
    band = []
    cmap = doppler_colormap()
    norm = Normalize(vmin=-2.5, vmax=2.5)
    rgb = cmap(norm(radvel))
    indnodata = np.where(validity == 0)
    for ich in range(3):
        channel = np.round(rgb[:, :, ich]*255).astype('uint8')
        channel[indnodata] = 0
        band.append({'array':channel, 'nodatavalue':0})
    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
