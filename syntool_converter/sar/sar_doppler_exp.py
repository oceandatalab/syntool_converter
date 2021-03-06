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
import os
import numpy as np
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
import numpy.ma as ma
from scipy.ndimage.filters import median_filter


def descalloping(radvel, sweepangle):
    """
    """
    # x = 'rvlSweepAngle'
    # corr = ...
    # -(20*x)*(abs(x) le .1)
    # -(8*(x+.15))*(x gt .1 and x le .54)
    # +(50*(x-.65))*(x gt .54)
    # -(8*(x-.15))*(x lt -.1 and x ge -.54)
    # +(50*(x+.65))*(x lt -.54)
    # corr = corr/106.
    # rvl_corrigé = rvl + corr
    corr = - 20 * sweepangle * (abs(sweepangle) <= .1) \
           - 8 * (sweepangle + .15) * ((sweepangle > .1) & (sweepangle <= .54)) \
           + 50 * (sweepangle - .65) * (sweepangle > .54) \
           - 8 * (sweepangle - .15) * ((sweepangle < -.1) & (sweepangle >= -.54)) \
           + 50 * (sweepangle + .65) * (sweepangle < -.54)
    corr /= 106
    return radvel + corr


def smooth(radvel, size=(3, 3)):
    """
    """
    return median_filter(radvel, size=size)


def sar_doppler_exp(infile, outdir, pngkml=False,
                    vmin=-2.5, vmax=2.5, vmin_pal=-2.5, vmax_pal=2.5):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    sardop = Dataset(infile)
    mission = sardop.MISSIONNAME
    if mission == 'S1A':
        sensor_name = 'Sentinel-1A'
        sensor_platform = 'Sentinel-1A'
        source_provider = 'ESA'
    else:
        raise Exception('S1A mission expected.')
    doptime = sardop.variables['rvlZeroDopplerTime'][:]
    start_time = datetime.strptime(''.join(list(doptime[0, 0, :])),
                                   '%Y-%m-%dT%H:%M:%S.%f')
    stop_time = datetime.strptime(''.join(list(doptime[-1, -1, :])),
                                  '%Y-%m-%dT%H:%M:%S.%f')
    heading = sardop.variables['rvlHeading'][:]
    if np.sin((90 - heading.mean()) * np.pi / 180) > 0:
        sensor_pass = 'Ascending'
    else:
        sensor_pass = 'Descending'
    # safe_name = os.path.basename(os.path.dirname(os.path.dirname(infile)))
    # sensor_mode = safe_name.split('_')[1]
    # if sensor_mode not in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'IW', 'EW']:
    #     raise Exception('S[1-6]/IW/EW modes expected.')
    # sensor_swath = os.path.basename(infile).split('-')[1].upper()
    # sensor_polarisation = sardop.read_global_attribute('polarisation')
    # datagroup = safe_name.replace('.SAFE', '')
    # pid = datagroup.split('_')[-1]
    # dataname = os.path.splitext(os.path.basename(infile))[0] + '-' + pid
    dataname = os.path.splitext(os.path.basename(infile))[0]
    sensor_mode = dataname.split('_')[1]
    sensor_swath = sensor_mode
    sensor_polarisation = sardop.POLARISATION
    radvel = sardop.variables['rvlRadVel'][:]
    sweepangle = sardop.variables['rvlSweepAngle'][:]
    radvel = descalloping(radvel, sweepangle)
    radvel = smooth(radvel)
    inc = sardop.variables['rvlIncidenceAngle'][:]
    radvel /= np.sin(np.deg2rad(inc))
    #landflag = sardop.variables['rvlLandFlag'][:]
    lon = sardop.variables['rvlLon'][:]
    lat = sardop.variables['rvlLat'][:]
    if sensor_pass == 'Ascending':
        radvel *= -1
    ngcps = np.ceil(np.array(lon.shape) / 10.) + 1
    pix = np.linspace(0, lon.shape[1] - 1, num=ngcps[1]).round().astype('int32')
    lin = np.linspace(0, lon.shape[0] - 1, num=ngcps[0]).round().astype('int32')
    pix2d, lin2d = np.meshgrid(pix, lin)
    gcplon = lon[lin2d, pix2d]
    gcplat = lat[lin2d, pix2d]
    gcppix = pix2d + 0.5
    gcplin = lin2d + 0.5
    gcphei = np.zeros(ngcps)
    if gcplon.min() < -135 and gcplon.max() > 135:
        gcplon[np.where(gcplon < 0)] += 360.

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SAR_doppler_exp'
    metadata['name'] = dataname
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = source_provider
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'radial horizontal velocities'
    metadata['type'] = 'remote sensing'
    metadata['sensor_type'] = 'SAR'
    metadata['sensor_name'] = sensor_name
    metadata['sensor_platform'] = sensor_platform
    metadata['sensor_mode'] = sensor_mode
    metadata['sensor_swath'] = sensor_swath
    metadata['sensor_polarisation'] = sensor_polarisation
    metadata['sensor_pass'] = sensor_pass
    # metadata['datagroup'] = datagroup
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                gcppix, gcplin)
    band = []
    #indndv = np.where(landflag != 0)
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(radvel, vmin, vmax, out=radvel)
    array = np.round((radvel - offset) / scale).astype('uint8')
    #array[indndv] = 255
    colortable = stfmt.format_colortable('doppler', vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'radial horizontal velocities', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    # Write projected png/kml
    if pngkml == True:
        print 'Write projected png/kml'
        stfmt.write_pngkml_proj(tifffile)
