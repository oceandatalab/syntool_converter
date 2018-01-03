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

from cerbere.mapper.safegeotifffile import SAFEGeoTiffFile
from sar.data.sarimage import SARImage
import numpy as np
from datetime import datetime
import os
import syntool_converter.utils.syntoolformat as stfmt
import re
import pickle as pkl
from scipy.stats import scoreatpercentile
import pyproj


def get_landmask(lon, lat, landmaskpath):
    """
    """
    from ceraux.landmask import Landmask
    from cerbere.datamodel.image import Image
    landmaskinstance = Landmask(filename=landmaskpath)
    feature = Image(longitudes=lon, latitudes=lat, times=datetime(1, 1, 1))
    landmaskfield = landmaskinstance.get_land_mask(feature)
    landmask = landmaskfield.get_values()
    # landmask should be a ndarray, use np.ma.getdata() just in case
    return np.ma.getdata(landmask)


def get_caltmp(fname, datet):
    """
    """
    cal = pkl.load(open(fname))
    deltat = abs(cal['time']-datet)
    caltmp = cal['correction_factor'][deltat.argmin()]
    return caltmp


# S1 resolutions (rng x azi) / pixel spacings (rng x azi)
# GRD
# WV GRDM 52 x 51 / 25 x 25
# SM GRDF 9 x 9 / 4 x 4
# SM GRDH 23 x 23 / 10 x 10
# SM GRDM 84 x 84 / 40 x 40
# IW GRDH 20 x 22 / 10 x 10
# IW GRDM 88 x 89 / 40 x 40
# EW GRDH 50 x 50 / 25 x 25
# EW GRDM 93 x 87 / 40 x 40
# SLC
# WV SLC 2.0 x 4.8 and 3.1 x 4.8 / 1.7 x 4.1 and 2.7 x 4.1
# SM SLC 1.7 x 4.3 to 3.6 x 4.9 / 1.5 x 3.6 to 3.1 x 4.1
# IW SLC 2.7 x 22 to 3.5 x 22 / 2.3 x 17.4 to 3 x 17.4
# EW SLC 7.9 x 42 to 14.4 x 43 / 9 x 34.7 to 12.5 x 34.7


def sar_roughness(infile, outdir, pngkml=False, contrast=None, vmin=None,
                  vmax=None, landmaskpath=None, write_netcdf=False,
                  gcp2height=0):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    sarmp = SAFEGeoTiffFile(infile)
    sarim = SARImage(sarmp)
    mission = sarim.get_info('mission')
    if mission == 'S1A':
        sensor_name = 'Sentinel-1A'
        sensor_platform = 'Sentinel-1A'
        source_provider = 'ESA'
    elif mission == 'S1B':
        sensor_name = 'Sentinel-1B'
        sensor_platform = 'Sentinel-1B'
        source_provider = 'ESA'
    else:
        raise Exception('Unknown mission')
    timefmt = '%Y-%m-%dT%H:%M:%S.%f'
    start_time = datetime.strptime(sarim.get_info('start_time'), timefmt)
    stop_time = datetime.strptime(sarim.get_info('stop_time'), timefmt)
    sensor_pass = sarim.get_info('pass')
    sensor_mode = sarim.get_info('mode')
    sensor_swath = sarim.get_info('swath')
    sensor_polarisation = sarim.get_info('polarisation')
    product = sarim.get_info('product')
    if product == 'GRD':
        spacing = [2, 2]
    elif product == 'SLC':
        if sensor_mode == 'WV':
            mspacing = (15, 15)
        elif re.match(r'^S[1-6]$', sensor_mode) != None:
            mspacing = (15, 15)
        elif sensor_mode == 'IW':
            raise Exception('sar_roughness for IW SLC ?')
        elif sensor_mode == 'EW':
            raise Exception('sar_roughness for EW SLC ?')
        else:
            raise Exception('Unkown S1 mode : {}'.format(sensor_mode))
        spacing = np.round(sarim.meters2pixels(mspacing))
    else:
        raise Exception('Unkown S1 product : {}'.format(product))
    mspacing = sarim.pixels2meters(spacing)
    datagroup = sarim.get_info('safe_name').replace('.SAFE', '')
    pid = datagroup.split('_')[-1]
    dataname = os.path.splitext(os.path.basename(infile))[0] + '-' + pid
    ssr = np.sqrt(sarim.get_data('roughness', spacing=spacing))
    ########## TMP calibration constant ##########
    # if sensor_mode == 'WV':
    #     caldir = '/home/cercache/project/mpc-sentinel1/analysis/s1_data_analysis/L1/WV/S1A_WV_SLC__1S/cal_cste'
    #     if sensor_polarisation == 'HH':
    #         if sensor_swath == 'WV1':
    #             caltmp = (55.80+56.91)/2.
    #             calname = 'cal_cste_hh_wv1.pkl'
    #         elif sensor_swath == 'WV2':
    #             caltmp = (40.65+40.32)/2.
    #             calname = 'cal_cste_hh_wv2.pkl'
    #     elif sensor_polarisation == 'VV':
    #         if sensor_swath == 'WV1':
    #             caltmp = 58.24
    #             calname = 'cal_cste_vv_wv1.pkl'
    #         elif sensor_swath == 'WV2':
    #             caltmp = 49.02
    #             calname = 'cal_cste_vv_wv2.pkl'
    #     calpath = os.path.join(caldir, calname)
    #     if os.path.exists(calpath) == True:
    #         caltmp = get_caltmp(calpath, start_time)
    # elif re.match(r'^S[1-6]$', sensor_mode) != None:
    #     if start_time < datetime(2014, 7, 16, 0, 0, 0):
    #         if sensor_mode == 'S6':
    #             raise Exception('S6 calibration missing')
    #         sm2cal = {'S1':58., 'S2':56., 'S3':52., 'S4':52., 'S5':49.}
    #     else:
    #         # from commissioning phase report
    #         sm2cal = {'S1':3., 'S2':5., 'S3':-1.5, 'S4':4., 'S5':1., 'S6':4.75}
    #     caltmp = sm2cal[sensor_mode]
    # elif sensor_mode == 'IW':
    #     if start_time < datetime(2014, 7, 16, 0, 0, 0):
    #         caltmp = 109.
    #     else:
    #         caltmp = 3. # from commissioning phase report
    # elif sensor_mode == 'EW':
    #     if start_time < datetime(2014, 7, 16, 0, 0, 0):
    #         caltmp = 94.
    #     else:
    #         caltmp = -1. # <- -2. # from commissioning phase report
    # else:
    #     raise Exception('Which tmp calibration constant for this mode ?')
    # print '--> caltmp=%f' % caltmp
    # ssr *= np.sqrt(10 ** (caltmp / 10.))
    ########## /TMP calibration constant ##########
    dim = ssr.shape
    # Set contrast
    if vmin == None or vmax == None:
        if contrast == None:
            if sensor_mode == 'WV':
                contrast = 'relative'
            else:
                contrast = 'sea'
        if contrast == 'relative':
            if sensor_mode == 'WV':
                noborder = [slice(int(dim[0]*.05), int(dim[0]*.95)),
                            slice(int(dim[1]*.05), int(dim[1]*.95))]
            else:
                noborder = [slice(int(dim[0]*.05), int(dim[0]*.95)),
                            slice(int(dim[1]*.1), int(dim[1]*.9))]
            values = ssr[noborder]
            if landmaskpath != None and os.path.exists(landmaskpath):
                lmspacing = np.round(sarim.meters2pixels(111.32/120*1000))
                lmspacing -= np.mod(lmspacing, spacing)
                lon = sarim.get_data('lon', spacing=lmspacing)
                lat = sarim.get_data('lat', spacing=lmspacing)
                lmdim = (lon.shape[0]+1, lon.shape[1]+1)
                landmask = np.ones(lmdim, dtype=bool)
                landmask[:-1, :-1] = get_landmask(lon, lat, landmaskpath)
                lmfac = lmspacing / spacing
                landmask = np.repeat(landmask, lmfac[0], axis=0)
                landmask = np.repeat(landmask, lmfac[1], axis=1)
                seaindex = np.where(landmask[noborder] == False)
                if seaindex[0].size >= ssr.size*0.01:
                    values = values[seaindex]
            if vmin == None:
                vmin = scoreatpercentile(values, 0.1)
            if vmax == None:
                vmax = scoreatpercentile(values, 99.9)
        elif contrast == 'sea':
            if sensor_polarisation in ['HH', 'VV']:
                if vmin == None:
                    vmin = 0.
                if vmax == None:
                    vmax = 2.
            else:
                if vmin == None:
                    vmin = 1.
                if vmax == None:
                    vmax = 3.
        elif contrast == 'ice':
            if sensor_polarisation in ['HH', 'VV']:
                if vmin == None:
                    vmin = 0.
                if vmax == None:
                    vmax = 3.5
            else:
                if vmin == None:
                    vmin = 1.
                if vmax == None:
                    vmax = 5.
        else:
            raise Exception('Unknown contrast name.')
    print '--> vmin=%f vmax=%f' % (vmin, vmax)
    ssr = ssr[::-1, :] # keep SAR orientation for geotiff
    geoloc = sarim.get_info('geolocation_grid')
    gcplin = (dim[0] * spacing[0] - 1 - geoloc['line'] + 0.5) / spacing[0]
    gcppix = (geoloc['pixel'] + 0.5) / spacing[1]
    gcplon = geoloc['longitude']
    gcplat = geoloc['latitude']
    gcphei = geoloc['height']
    if gcp2height is not None:
        geod = pyproj.Geod(ellps='WGS84')
        gcpforw, gcpback, _ = geod.inv(gcplon[:, :-1], gcplat[:, :-1],
                                       gcplon[:, 1:], gcplat[:, 1:])
        gcpforw = np.hstack((gcpforw, gcpforw[:, [-1]]))
        gcpback = np.hstack((gcpback[:, [0]], gcpback))
        gcpinc = geoloc['incidence_angle']
        mvdist = (gcp2height - gcphei) / np.tan(np.deg2rad(gcpinc))
        mvforw = gcpforw
        indneg = np.where(mvdist < 0)
        mvdist[indneg] = -mvdist[indneg]
        mvforw[indneg] = gcpback[indneg]
        _gcplon, _gcplat, _ = geod.fwd(gcplon, gcplat, mvforw, mvdist)
        gcplon = _gcplon
        gcplat = _gcplat
        gcphei.fill(gcp2height)
    if gcplon.min() < -135 and gcplon.max() > 135:
        gcplon[np.where(gcplon < 0)] += 360.
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    if sensor_polarisation in ['HH', 'VV']:
        metadata['product_name'] = 'SAR_roughness'
    else:
        metadata['product_name'] = 'SAR_roughness_crosspol'
    metadata['name'] = dataname
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = source_provider
    metadata['processing_center'] = 'OceanDataLab'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea surface roughness'
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
    geolocation['projection'] = sarim._mapper._handler.GetGCPProjection()
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                gcppix, gcplin)
    band = []
    scale = (vmax-vmin)/254.
    offset = vmin
    indzero = np.where(ssr == 0)
    array = np.clip(np.round((ssr-offset)/scale), 0, 254).astype('uint8')
    array[indzero] = 255
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea surface roughness', 'unittype':'',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax]})
    # Write
    if write_netcdf == False:
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
        # Write projected png/kml
        if pngkml == True:
            print 'Write projected png/kml'
            stfmt.write_pngkml_proj(tifffile)
    elif write_netcdf == True:
        print 'Write netcdf'
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        band[0]['name'] = 'sea_surface_roughness'
        band[0]['long_name'] = 'sea surface roughness'
        band[0]['unittype'] = '1'
        metadata['spatial_resolution'] = mspacing.min()
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'swath',
                           ngcps=gcplon.shape)
