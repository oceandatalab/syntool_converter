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


def sar_wind(infile, outdir, pngkml=False, valid_percent_min=1.,
             vmin=0., vmax=25.4, vmin_pal=0., vmax_pal=50*0.514):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    sarwind = SAFEOCNNCFile(infile, product='WIND')
    mission = sarwind.read_global_attribute('missionName')
    if mission == 'S1A':
        sensor_name = 'Sentinel-1A'
        sensor_platform = 'Sentinel-1A'
        source_provider = 'ESA'
    elif mission == 'S1B':
        sensor_name = 'Sentinel-1B'
        sensor_platform = 'Sentinel-1B'
        source_provider = 'ESA'
    else:
        raise Exception('S1A/S1B missions expected.')
    start_time = sarwind.get_start_time()
    stop_time = sarwind.get_end_time()
    heading = sarwind.read_values('owiHeading')
    if np.sin((90 - heading[0, 0]) * np.pi / 180) > 0:
        sensor_pass = 'Ascending'
    else:
        sensor_pass = 'Descending'
    safe_name = os.path.basename(os.path.dirname(os.path.dirname(infile)))
    sensor_mode = safe_name.split('_')[1]
    if sensor_mode not in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'IW', 'EW']:
        raise Exception('S[1-6]/IW/EW modes expected.')
    sensor_swath = os.path.basename(infile).split('-')[1].upper()
    sensor_polarisation = sarwind.read_global_attribute('polarisation')
    datagroup = safe_name.replace('.SAFE', '')
    pid = datagroup.split('_')[-1]
    dataname = os.path.splitext(os.path.basename(infile))[0] + '-' + pid
    windspeed = sarwind.read_values('owiWindSpeed')
    if windspeed.shape == (1, 1):
        raise Exception('owiRaSize and owiAzSize equals 1 !')
    winddirection = sarwind.read_values('owiWindDirection')
    landflag = sarwind.read_values('owiLandFlag')
    inversionquality = sarwind.read_values('owiInversionQuality')
    windquality = sarwind.read_values('owiWindQuality')
    #pbright = sarwind.read_values('owiPBright')
    lon = sarwind.read_values('lon')
    lat = sarwind.read_values('lat')
    if np.ma.is_masked(lon) or np.ma.is_masked(lat):
        raise Exception('Some lon and/or lat is masked.')
    if np.all(lon == 0) or np.all(lat == 0):
        raise Exception('All lon and/or lat set to 0.')
    ngcps = np.ceil(np.array(lon.shape) / 10.).astype('int') + 1
    pix = np.linspace(0, lon.shape[1] - 1, num=ngcps[1]).round().astype('int32')
    lin = np.linspace(0, lon.shape[0] - 1, num=ngcps[0]).round().astype('int32')
    pix2d, lin2d = np.meshgrid(pix, lin)
    gcplon = lon[lin2d, pix2d]
    gcplat = lat[lin2d, pix2d]
    gcppix = pix2d + 0.5
    gcplin = lin2d + 0.5
    gcphei = np.zeros(ngcps)
    ## Make sure lon are continuous (no jump because of IDL crossing)
    ## (if IDL crossing, by convention we make lon to be around 180deg)
    # if gcplon.min() < -135 and gcplon.max() > 135:
    #     gcplon[np.where(gcplon < 0)] += 360.
    gcplonmid = gcplon[ngcps[0] / 2, ngcps[1] / 2]
    gcplon = np.mod(gcplon - (gcplonmid - 180.), 360.) + gcplonmid - 180.
    gcplonmin = gcplon.min()
    gcplon = gcplon - np.floor((gcplonmin + 180.) / 360.) * 360.

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SAR_wind'
    metadata['name'] = dataname
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = source_provider
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = ['wind speed', 'wind direction']
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
    #mask = landflag != 0
    mask = (landflag != 0) | \
        ((windspeed == 0) & (winddirection == 180)) | \
        ((windspeed == 0) & (windquality == 3)) | \
        ((windspeed == 0) & (inversionquality == 2))
    mask = np.ma.getdata(mask) # we don't want to sum on a masked mask
    valid_percent = np.sum(~mask) / float(mask.size) * 100
    if valid_percent <= valid_percent_min:
        raise Exception('Not enough valid data: {:0.3f}%'.format(valid_percent))
    # if np.all(mask):
    #     raise Exception('Data is fully masked !')
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(windspeed, vmin, vmax, out=windspeed)
    array = np.round((windspeed - offset) / scale).astype('uint8')
    array[mask] = 255
    colortable = stfmt.format_colortable('noaa_wind', vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'wind speed', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    winddirection = np.mod(90. - winddirection + 180., 360.)
    array = np.round(winddirection/360.*254.).astype('uint8')
    array[mask] = 255
    band.append({'array':array, 'scale':360./254., 'offset':0.,
                 'description':'wind direction', 'unittype':'deg',
                 'nodatavalue':255, 'parameter_range':[0, 360.]})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    # Write projected png/kml
    if pngkml == True:
        print 'Write projected png/kml'
        stfmt.write_pngkml_proj(tifffile)
