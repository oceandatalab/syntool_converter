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

from cerbere.mapper.safeocnncfile import SAFEOCNNCFile
import os
import numpy as np
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
import numpy.ma as ma


def sar_doppler(infile, outdir, pngkml=False,
                vmin=-2.5, vmax=2.5, vmin_pal=-2.5, vmax_pal=2.5):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    sardop = SAFEOCNNCFile(infile, product='DOPPLER')
    mission = sardop.read_global_attribute('missionName')
    if mission == 'S1A':
        sensor_name = 'Sentinel-1A'
        sensor_platform = 'Sentinel-1A'
        source_provider = 'ESA'
    else:
        raise Exception('S1A mission expected.')
    start_time = sardop.get_start_time()
    stop_time = sardop.get_end_time()
    heading = sardop.read_values('rvlHeading')
    if np.sin((90 - heading.mean()) * np.pi / 180) > 0:
        sensor_pass = 'Ascending'
    else:
        sensor_pass = 'Descending'
    safe_name = os.path.basename(os.path.dirname(os.path.dirname(infile)))
    sensor_mode = safe_name.split('_')[1]
    if sensor_mode not in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'IW', 'EW']:
        raise Exception('S[1-6]/IW/EW modes expected.')
    sensor_swath = os.path.basename(infile).split('-')[1].upper()
    sensor_polarisation = sardop.read_global_attribute('polarisation')
    datagroup = safe_name.replace('.SAFE', '')
    pid = datagroup.split('_')[-1]
    dataname = os.path.splitext(os.path.basename(infile))[0] + '-' + pid
    if 'rvlSwath' in sardop.get_dimensions():
        nswath = sardop.get_dimsize('rvlSwath')
    else:
        nswath = 1
    for iswath in range(nswath):

        if nswath == 1:
            radvel = sardop.read_values('rvlRadVel')
            landflag = sardop.read_values('rvlLandFlag')
            lon = sardop.read_values('lon')
            lat = sardop.read_values('lat')
            name = dataname
        else:
            radvel = sardop.read_values('rvlRadVel')[:, :, iswath]
            landflag = sardop.read_values('rvlLandFlag')[:, :, iswath]
            lon = sardop.read_values('lon')[:, :, iswath]
            lat = sardop.read_values('lat')[:, :, iswath]
            valid = np.where((ma.getmaskarray(lon) == False) & \
                             (ma.getmaskarray(lat) == False))
            slices = [slice(valid[0].min(), valid[0].max() + 1),
                      slice(valid[1].min(), valid[1].max() + 1)]
            radvel = radvel[slices]
            landflag = landflag[slices]
            lon = lon[slices]
            lat = lat[slices]
            name = dataname + '-' + str(iswath+1)

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
        metadata['product_name'] = 'SAR_doppler'
        metadata['name'] = name
        metadata['datetime'] = dtime
        metadata['time_range'] = time_range
        metadata['source_URI'] = infile
        metadata['source_provider'] = source_provider
        metadata['processing_center'] = '' #'OceanDataLab'
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
        metadata['datagroup'] = datagroup
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
