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
import numpy as np
import numpy.ma as ma
from scipy.ndimage.morphology import binary_opening
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import syntool_converter.utils.resample as resample
from syntool_converter.utils.sst import denoise_sst
import os
import pyproj
import sys

#import geometry.point


# swath size : ~4044.00km x ~3114.07km
# 100km ngcps=(41, 32)
# viirs_sst('/mnt/data/viirs/20150303113959-OSPO-L2P_GHRSST-SST1m-VIIRS_NPP-v02.0-fv01.0.nc',
#           '/mnt/data/syntool_inputs/synergy/')

def viirs_sst(infile, outdir, vmin=None, vmax=None, contrast='relative',
              ngcps=(41, 32), denoise_kernel='boxcar', denoise_width=27,
              open_iterations=1, nprocs=1,
              pngkml=False, write_netcdf=False, file_range=None):
    """
    """
    dic = None

    # Check file containing ranges
    if file_range is not None:
        if not os.path.isfile(file_range):
            raise Exception('file_range {} not found'.format(file_range))
        # Read a txt file which contains three columns: yearday,vmin,vmax
        with open(file_range, 'r') as f:
            dic = {}
            for line in f:
                (fdoy, fmin, fmax) = line.split(',')
                dic[int(fdoy)] = (float(fmin), float(fmax))

    if contrast == 'med':
        listbox = [[-6., 35., 2.75, 42.48],
                   [2.74, 30, 42.2, 47.00]]
    elif contrast == 'cwe':
        listbox = [[-23., 35.2, -5.5, 42.88],
                   [-23., 42.8, 2.20, 51.]]
    elif contrast == 'nwe':
        listbox = [[-23., 50.8, 32.7, 68.]]
    elif contrast == 'gom':
        listbox = [[-98., 18.0, -80.5, 30.5]]
    elif contrast == 'agulhas':
        listbox = [[10.8437, -45.7404, 39.9799, -25.3019]]
    elif contrast == 'gs':
        listbox = [[-81.52, 20, -30, 45]]
    else:
        listbox = None
    # Read/Process data
    print 'Read/Process data'
    dataset = Dataset(infile)
    start_time = datetime.strptime(dataset.start_time, '%Y%m%dT%H%M%SZ')
    print start_time.day
    print start_time.month
    stop_time = datetime.strptime(dataset.stop_time, '%Y%m%dT%H%M%SZ')
    lon = dataset.variables['lon'][:, :]
    lat = dataset.variables['lat'][:, :]
    sst = np.ma.array(dataset.variables['sea_surface_temperature'][0, :, :])
    _bt11= dataset.variables['brightness_temperature_11um'][0, :, :]
    bt11 = np.ma.array(_bt11)
    quality_level = np.ma.array(dataset.variables['quality_level'][0, :, :])
    '''
    if file_shape is not None:
        with open(file_shape, 'r') as fshape:
            shape = shapely.wkt.load(fshape)
        box = shape.bounds
        index_in = np.where((lon >= box[0]) & (lat >= box[1])
                            & (lon <= box[2]) & (lat <= box[3]))
        index_out = np.where((lon < box[0]) | (lat < box[1])
                             | (lon > box[2]) | (lat > box[3]))
        sst[index_out] = np.nan
        print(np.shape(index_in))
        sys.exit(1)
        for i, j in zip(index_in[0], index_in[1]):
            p = Point(lon[i, j], lat[i, j])
            if p.within(shape) is False:
                sst[i, j] = np.nan
    '''
    if listbox is not None:
        mask_box = np.zeros(np.shape(sst))
        for i in range(np.shape(listbox)[0]):
            index_in = np.where((lon >= listbox[i][0]) & (lat >= listbox[i][1])
                             & (lon <= listbox[i][2]) & (lat <= listbox[i][3]))
            mask_box[index_in] = 1
        mask = ma.getmaskarray(sst) | ma.getmaskarray(bt11) | \
               (quality_level.data < 4) | (mask_box == 0)
    else:
        mask = ma.getmaskarray(sst) | ma.getmaskarray(bt11) | \
               (quality_level.data < 4)
    if mask.all():
        print 'No data'
        sys.exit(0)
    # GCPs for resampling and geotiff georeference
    scansize = 16
    dtime0 = datetime.utcnow()
    gcps = resample.get_gcps_from_bowtie(lon, lat, scansize, ngcps=ngcps)
    dtime = datetime.utcnow() - dtime0
    print 'Get GCPs from bowtie swath : {}'.format(dtime)
    gcplon, gcplat, gcpnpixel, gcpnline = gcps
    rspysize = lon.shape[0]
    geod = pyproj.Geod(ellps='WGS84')
    mid = abs(gcpnline[:, 0] - 0.5).argmin()
    xdists = geod.inv(gcplon[mid, :-1], gcplat[mid, :-1],
                      gcplon[mid, 1:], gcplat[mid, 1:])[2]
    xdist = np.sum(xdists) / abs(gcpnpixel[mid, -1] - gcpnpixel[mid, 0])
    rspxsize = np.round(xdist / 750.).astype('int') + 1
    gcpline = gcpnline * rspysize
    gcppixel = gcpnpixel * rspxsize

    # Resample with LinearNDInterpolator in output space
    dtime0 = datetime.utcnow()
    pix, lin = resample.get_points_from_gcps(gcplon, gcplat, gcppixel,
                                             gcpline, rspxsize, rspysize,
                                             1, lon, lat, nprocs=nprocs) - 0.5
    dtime = datetime.utcnow() - dtime0
    print 'Get input coordinates in new grid : {}'.format(dtime)
    # Test input grid in output space
    # import matplotlib.pyplot as plt
    # for iscan in range(lon.shape[0] / scansize):
    #     pixscan = pix[iscan * scansize: (iscan+1) * scansize, :]
    #     linscan = lin[iscan * scansize: (iscan+1) * scansize, :]
    #     # maskscan = mask[iscan * scansize: (iscan+1) * scansize, :]
    #     # pixscan = pixscan[~maskscan]
    #     # linscan = linscan[~maskscan]
    #     plt.plot(pixscan.flatten(), linscan.flatten(), '+')
    # plt.show()
    # import pdb ; pdb.set_trace()
    # \Test input grid in output space
    dtime0 = datetime.utcnow()
    sst.data[mask] = np.nan
    bt11.data[mask] = np.nan
    val = np.dstack((sst.data, bt11.data))
    rspval = resample.resample_bowtie_linear(pix, lin, val, scansize,
                                             rspxsize, rspysize, show=False)
    rspsst = rspval[:, :, 0]
    rspbt11 = rspval[:, :, 1]
    rspmask = ma.getmaskarray(rspsst) | ma.getmaskarray(rspbt11)
    dtime = datetime.utcnow() - dtime0
    print 'Interpolate in new grid : {}'.format(dtime)

    # Denoise sst and open mask
    rspsst.mask = rspmask
    rspbt11.mask = rspmask
    finalsst = denoise_sst(rspsst, rspbt11, kernel=denoise_kernel,
                           width=denoise_width, show=False)
    finalmask = ~binary_opening(~rspmask, structure=np.ones((3, 3)),
                                iterations=open_iterations)
    finalsst.mask = finalmask

    # Contrast
    if vmin == None:
        if contrast == 'relative':
            vmin = np.percentile(finalsst.compressed(), 0.5)
        #elif contrast == 'agulhas':
        #    dayofyear = float(start_time.timetuple().tm_yday)
        #    vmin = 273.15 + 2. * np.cos((dayofyear - 45.) * 2. * np.pi / 365.) + 20. - 9.
        #    #par = [277.94999694824219, 42, 2.5500030517578125, -219]
        #    par = [278.09999084472656, 0.62831853071795862,
        #           2.4000091552734375, 0.1570796326794896]
        #    vmin = par[0] + par[2] * np.cos(par[3] * dayofyear - par[1])
        #if a specific txt file is provided for the range
        elif dic is not None:
            dayofyear = float(start_time.timetuple().tm_yday)
            extrema = dic.get(dayofyear, dic[min(dic.keys(),
                           key=lambda k:abs(k - dayofyear))])
            vmin = extrema[0]
        else:
            raise Exception('Unknown contrast : {}'.format(contrast))
    if vmax == None:
        if contrast == 'relative':
            vmax = np.percentile(finalsst.compressed(), 99.5)
        #elif contrast == 'agulhas':
        #    dayofyear = float(start_time.timetuple().tm_yday)
        #    vmax = 273.15 + 2. * np.cos((dayofyear - 45.) * 2. * np.pi / 365.) + 20. + 4.
        #    #par = [300.59999084472656, 21, 2.8499908447265625, -191]
        #    par = [300.59999084472656, 0.29919930034188508,
        #           2.8499908447265625, 0.14959965017094254]
        #    vmax = par[0] + par[2] * np.cos(par[3] * dayofyear - par[1])
        #if a specific text file is provided for the range
        elif dic is not None:
            dayofyear = float(start_time.timetuple().tm_yday)
            extrema = dic.get(dayofyear, dic[min(dic.keys(),
                           key=lambda k:abs(k - dayofyear))])
            vmax = extrema[1]
        else:
            raise Exception('Unknown contrast : {}'.format(contrast))

    # Flip (geotiff in "swath sense")
    finalsst = finalsst[::-1, ::-1]
    gcppixel = rspxsize - gcppixel
    gcpline = rspysize - gcpline

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SST_VIIRS_denoised'
    if contrast == 'relative':
        metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    else:
        metadata['name'] = '{}_{}'.format(os.path.splitext(os.path.basename(infile))[0], contrast)
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'NOAA'
    metadata['processing_center'] = 'OceanDataLab'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea surface temperature'
    metadata['type'] = 'remote sensing'
    metadata['sensor_type'] = 'radiometer'
    metadata['sensor_name'] = 'VIIRS'
    metadata['sensor_platform'] = 'Suomi-NPP'
    #metadata['sensor_pass'] =
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    gcpheight = np.zeros(gcppixel.shape)
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcpheight,
                                                gcppixel, gcpline)
    band = []
    indndv = np.where(ma.getmaskarray(finalsst) == True)
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(finalsst.data, vmin, vmax, out=finalsst.data)
    array = np.round((finalsst.data - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('cerbere_medspiration',
                                         vmax=vmax, vmax_pal=vmax,
                                         vmin=vmin, vmin_pal=vmin)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'sea surface temperature', 'unittype':'K',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    if write_netcdf == False:
        # Write geotiff
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
        band[0]['name'] = 'denoised_sst'
        band[0]['long_name'] = 'denoised sea surface temperature'
        band[0]['standard_name'] = 'sea_surface_temperature'
        # ymid = abs(gcpline[:, 0] - rspysize / 2.).argmin()
        # xdists = geod.inv(gcplon[ymid, :-1], gcplat[ymid, :-1],
        #                   gcplon[ymid, 1:], gcplat[ymid, 1:])[2] / \
        #                   np.abs(gcppixel[ymid, 1:] - gcppixel[ymid, :-1])
        # xmid = abs(gcppixel[0, :] - rspxsize / 2.).argmin()
        # ydists = geod.inv(gcplon[:-1, xmid], gcplat[:-1, xmid],
        #                   gcplon[1:, xmid], gcplat[1:, xmid])[2] / \
        #                   np.abs(gcpline[1:, xmid] - gcpline[:-1, xmid])
        # print xdists.min(), xdists.max(), xdists.mean()
        # # e.g. 749.905437495 749.905892002 749.905827652
        # print ydists.min(), ydists.max(), ydists.mean()
        # # e.g. 737.638084996 741.195663083 739.157662785
        metadata['spatial_resolution'] = 750.
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'swath',
                           ngcps=gcplon.shape)

