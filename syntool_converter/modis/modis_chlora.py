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

import numpy as np
import numpy.ma as ma
from scipy.ndimage.morphology import binary_opening
from datetime import datetime
import syntool_converter.utils.modis as modis
import syntool_converter.utils.syntoolformat as stfmt
import syntool_converter.utils.resample as resample
import os
import pyproj
import re
import math
import sys

# modis_chlora('A2011338122500', '/mnt/data/syntool_inputs/synergy/',
#              download_dir='/mnt/data/modis')

def modis_chlora(infileid, outdir, download_dir='/tmp',
                 vmin=None, vmax=None, contrast='relative',
                 ngcps=(21, 25), open_iterations=1,
                 nprocs=1, pngkml=False, write_netcdf=False):
    """
    """
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
    # Search/Download data
    print 'Search/Download data'
    if re.match(r'^[AT][0-9]{13}$', infileid) is None:
        raise Exception('Input for modis_chlora is an ID '
                        '(e.g. A2011338122500 or T2014143234500)')
    platform = infileid[0]
    date = datetime.strptime(infileid[1:], '%Y%j%H%M%S')
    modisocid = {'A':'MODISAL2OC', 'T':'MODISTL2OC'}[platform]
    modisocfname = modis.search_and_download(modisocid, date, download_dir)
    modis03id = {'A':'MYD03', 'T':'MOD03'}[platform]
    modis03fname = modis.search_and_download(modis03id, date, download_dir)

    # Read/Process data
    print 'Read/Process data'
    # Read from OC file
    modisocfile = modis.MODISL2File(modisocfname)
    # lon = modisocfile.read_lon()
    # lat = modisocfile.read_lat()
    chlora = modisocfile.read_chlora()
    attrs = modisocfile.read_attributes()
    modisocfile.close()
    # Read from geolocation file
    modis03file = modis.MODIS03File(modis03fname)
    lon = modis03file.read_lon()
    lat = modis03file.read_lat()
    modis03file.close()
    if listbox is not None:
        mask_box = np.zeros(np.shape(chlora.data))
        for i in range(np.shape(listbox)[0]):
            index_in = np.where((lon >= listbox[i][0]) & (lat >= listbox[i][1])
                             & (lon <= listbox[i][2]) & (lat <= listbox[i][3]))
            mask_box[index_in] = 1
        mask = (mask_box == 0) | ma.getmaskarray(chlora)
    else:
        mask = ma.getmaskarray(chlora)
    if mask.all():
        print 'No data'
        sys.exit(0)
    # GCPs for resampling and geotiff georeference
    scansize = 10
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
    rspxsize = np.round(xdist / 1000.).astype('int') + 1
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
    chlora.data[mask] = np.nan
    rspchlora = resample.resample_bowtie_linear(pix, lin, chlora.data, scansize,
                                                rspxsize, rspysize, show=False)
    rspmask = ma.getmaskarray(rspchlora)
    dtime = datetime.utcnow() - dtime0
    print 'Interpolate in new grid : {}'.format(dtime)

    # Take log and open mask
    finalchlora = ma.log(rspchlora)
    finalmask = ~binary_opening(~rspmask, structure=np.ones((3, 3)),
                                iterations=open_iterations)
    finalchlora.mask = finalmask

    # Contrast
    if vmin == None:
        if contrast == 'relative':
            vmin = np.percentile(finalchlora.compressed(), 0.5)
        elif contrast == 'agulhas':
            dayofyear = float(attrs['start_time'].timetuple().tm_yday)
            vmin = -0.5 * np.cos((dayofyear - 45.) * 2. * np.pi / 365.) - 3.
        elif contrast == 'med' or contrast == 'nwe' or contrast == 'cwe':
            vmin = np.percentile(finalchlora.compressed(), 2)
        else:
            raise Exception('Unknown contrast : {}'.format(contrast))
    else:
        if vmin != 0:
            vmin = math.log(vmin)
        else:
            vmin = np.percentile(finalchlora.compressed(), 0.5)
    if vmax == None:
        if contrast == 'relative':
            vmax = np.percentile(finalchlora.compressed(), 99.5)
        elif contrast == 'agulhas':
            dayofyear = float(attrs['start_time'].timetuple().tm_yday)
            vmax = 0.5 * np.cos((dayofyear - 45.) * 2. * np.pi / 365.) + 3.
        elif contrast == 'med':
            vmax = np.percentile(finalchlora.compressed(), 98)
        elif contrast == 'nwe':
            vmax = np.percentile(finalchlora.compressed(), 98)
        elif contrast == 'cwe':
            vmax = np.percentile(finalchlora.compressed(), 98)
        else:
            raise Exception('Unknown contrast : {}'.format(contrast))
    else:
        if vmax != 0:
            vmax = math.log(vmax)
        else:
            vmax = np.percentile(finalchlora.compressed(), 98)

    # Flip (geotiff in "swath sense")
    finalchlora = finalchlora[::-1, ::-1]
    gcppixel = rspxsize - gcppixel
    gcpline = rspysize - gcpline

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(attrs['start_time'],
                                                      attrs['stop_time'],
                                                      units='ms')
    metadata['product_name'] = 'Chlorophyll_a_concentration_MODIS'
    if contrast == 'relative':
        metadata['name'] = os.path.splitext(os.path.basename(modisocfname))[0]
    else:
        metadata['name'] = '{}_{}'.format(os.path.splitext(os.path.basename(modisocfname))[0], contrast)
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = [modisocfname, modis03fname]
    metadata['source_provider'] = 'NASA'
    metadata['processing_center'] = 'OceanDataLab'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'chlorophyll a concentration'
    metadata['type'] = 'remote sensing'
    metadata['sensor_type'] = 'radiometer'
    metadata['sensor_name'] = 'MODIS'
    metadata['sensor_platform'] = attrs['platform']
    metadata['sensor_pass'] = attrs['pass']
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    gcpheight = np.zeros(gcppixel.shape)
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcpheight,
                                                gcppixel, gcpline)
    band = []
    indndv = np.where(ma.getmaskarray(finalchlora) == True)
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(finalchlora.data, vmin, vmax, out=finalchlora.data)
    array = np.round((finalchlora.data - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('chla_jet',
                                         vmax=vmax, vmax_pal=vmax,
                                         vmin=vmin, vmin_pal=vmin)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'chlorophyll a concentration',
                 'unittype':'log(mg/m3)',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
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
        band[0]['name'] = 'chlor_a'
        band[0]['long_name'] = 'Chlorophyll Concentration, OCI Algorithm'
        band[0]['standard_name'] = 'mass_concentration_chlorophyll_concentration_in_sea_water'
        band[0]['unittype'] = 'mg m^-3 (log)'
        # ymid = abs(gcpline[:, 0] - rspysize / 2.).argmin()
        # xdists = geod.inv(gcplon[ymid, :-1], gcplat[ymid, :-1],
        #                   gcplon[ymid, 1:], gcplat[ymid, 1:])[2] / \
        #                   np.abs(gcppixel[ymid, 1:] - gcppixel[ymid, :-1])
        # xmid = abs(gcppixel[0, :] - rspxsize / 2.).argmin()
        # ydists = geod.inv(gcplon[:-1, xmid], gcplat[:-1, xmid],
        #                   gcplon[1:, xmid], gcplat[1:, xmid])[2] / \
        #                   np.abs(gcpline[1:, xmid] - gcpline[:-1, xmid])
        # print xdists.min(), xdists.max(), xdists.mean()
        # # e.g. 999.569936185 999.569945796 999.569941619
        # print ydists.min(), ydists.max(), ydists.mean()
        # # e.g. 1005.25062313 1007.89500181 1006.72031127
        metadata['spatial_resolution'] = 1000.
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'swath',
                           ngcps=gcplon.shape)
