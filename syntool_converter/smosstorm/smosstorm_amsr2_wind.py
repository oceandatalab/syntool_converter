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
import numpy as np
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os


# def check_cont_lon():
#     """
#     """
#     import glob
#     paths = glob.glob('/mnt/data/smosstorm/data/smosstorm/l2/amsr2/2014/10/*/*.nc')
#     paths = ['/tmp/SOLAB_AMSR2_L2_NN_20160324_064324_20160324_073255_006_A_v1.nc']
#     for path in paths:
#         dataset = Dataset(path)
#         dataset.variables['longitude'].set_auto_mask(False)
#         lon = dataset.variables['longitude'][:]
#         if lon.min() < -180. or lon.max() > 180.:
#             raise Exception('Unexpected lon outside [-180., 180.]')
#         lon = lon.transpose()
#         lon_range = lon.max() - lon.min()
#         nrow, ncol = lon.shape
#         _nrow = np.ceil(nrow / 20.).astype('int')
#         _irow = np.linspace(0, nrow, num=_nrow, endpoint=False).round().astype('int')
#         _lon = lon.copy()
#         for irow in _irow:
#             lon0 = lon[irow, ncol / 2] - 180
#             lon[irow:, :] = np.mod(lon[irow:, :] - lon0, 360.) + lon0
#         lon_range2 = lon.max() - lon.min()
#         print lon_range, lon_range2, os.path.basename(path)
#         import matplotlib.pyplot as plt
#         lat = dataset.variables['latitude'][:]
#         lat = lat.transpose()
#         plt.scatter(_lon.flatten(), lat.flatten(), c='r', edgecolors='face')
#         plt.figure()
#         plt.scatter(lon.flatten(), lat.flatten(), c='b', edgecolors='face')
#         plt.show()


def smosstorm_amsr2_wind(infile, outdir, valid_max_lat=85.,
                         vmin=0., vmax=50.8, vmin_pal=0., vmax_pal=50.):
    """
    """
    # Read/Process data
    ## Get wind/lon/lat/time
    dataset = Dataset(infile)
    wind = dataset.variables['wind_speed_lf'][:, :]
    lat = dataset.variables['latitude'][:]
    if np.ma.is_masked(lat):
        raise Exception('Unexpected masked lat')
    ### For lon: _FillValue = -9999 with scale_factor = 0.01 makes -99.99 to be masked !
    ### -> we assume -99.99 is actually a correct longitude
    ### -> auto masking is disabled and we just check that lon are in [-180, 180]
    dataset.variables['longitude'].set_auto_mask(False)
    lon = dataset.variables['longitude'][:]
    if lon.min() < -180. or lon.max() > 180.:
        raise Exception('Unexpected lon outside [-180., 180.]')
    time = dataset.variables['time'][:]
    time_units = dataset.variables['time'].units
    ## Transpose: (col,row) -> (row,col)
    wind = wind.transpose()
    lat = lat.transpose()
    lon = lon.transpose()
    ## Keep only the minimum row slice with valid data
    ## (ie we don't keep first/last invalid data in along-track dimension)
    valid = ~np.ma.getmaskarray(wind)
    if valid_max_lat is not None:
        # We don't want to take into account points where lat > valid_max_lat
        valid[np.abs(lat) > valid_max_lat] = False
    row_valid = np.where(np.any(valid, axis=1))
    if row_valid[0].size == 0:
        raise Exception('Data is fully masked.')
    row_valid_min, row_valid_max = row_valid[0].min(), row_valid[0].max()
    if row_valid_min == row_valid_max: # For having at least 2 GCPs in row dim
        if row_valid_min == 0:
            row_valid_max += 1
        else:
            row_valid_min -= 1
    row_slice = slice(row_valid_min, row_valid_max + 1)
    #row_slice = slice(None) # for test
    rowcol_slice = [row_slice, slice(None)]
    wind = wind[rowcol_slice]
    lat = lat[rowcol_slice]
    lon = lon[rowcol_slice]
    time = time[row_slice]
    ## Make lon continuous
    nrow, ncol = lon.shape
    _nrow = np.ceil(nrow / 20.).astype('int')
    _irow = np.linspace(0, nrow, num=_nrow, endpoint=False).round().astype('int')
    for irow in _irow:
        lon0 = lon[irow, ncol / 2] - 180
        lon[irow:, :] = np.mod(lon[irow:, :] - lon0, 360.) + lon0
    lon = lon - np.floor((lon.min() + 180.) / 360.) * 360.
    lon_range = lon.max() - lon.min()
    #print lon.min(), lon.max(), lon_range
    if lon_range > 310.:
        raise Exception('Unexpected lon range : {}'.format(lon_range))
    ## Make GCPs
    dgcp = 20. # spacing is about 10km
    ngcps = np.ceil(np.array(lon.shape) / dgcp) + 1.
    lin = np.linspace(0, lon.shape[0] - 1, num=ngcps[0]).round().astype('int32')
    pix = np.linspace(0, lon.shape[1] - 1, num=ngcps[1]).round().astype('int32')
    pix2d, lin2d = np.meshgrid(pix, lin)
    pix2d = pix2d.flatten()
    lin2d = lin2d.flatten()
    ### Commented code : we assume now that lon are correct (_FillValue is not)
    # invalidgcp = np.where((ma.getmaskarray(lon[lin2d, pix2d]) == True) | \
    #                       (ma.getmaskarray(lat[lin2d, pix2d]) == True))
    # if invalidgcp[0].size != 0:
    #     validlonlat = np.where((ma.getmaskarray(lon) == False) & \
    #                            (ma.getmaskarray(lat) == False))
    #     for x in invalidgcp[0]:
    #         dists = np.sqrt((lin2d[x] - validlonlat[0]) ** 2 + \
    #                         (pix2d[x] - validlonlat[1]) ** 2)
    #         argmin = dists.argmin()
    #         if dists[argmin] <= dgcp / 2.:
    #             lin2d[x] = validlonlat[0][argmin]
    #             pix2d[x] = validlonlat[1][argmin]
    #         else:
    #             raise Exception('Invalid lon/lat for GCPs.')
    gcplon = lon[lin2d, pix2d]
    gcplat = lat[lin2d, pix2d]
    gcppix = pix2d + 0.5
    gcplin = lin2d + 0.5
    gcphei = np.zeros(gcplon.shape)
    ## Dataset name and time range
    ### Commented code : used with old smosstorm files (storm colocated)
    # head, tail = os.path.split(infile)
    # head, tail = os.path.split(head)
    # if tail != 'data':
    #     raise Exception('Input file parent directories are not '
    #                     '<storm_name>/amsr2/data !')
    # head, tail = os.path.split(head)
    # if tail != 'amsr2':
    #     raise Exception('Input file parent directories are not '
    #                     '<storm_name>/amsr2/data !')
    # head, storm_name = os.path.split(head)
    # dataset_name = os.path.splitext(os.path.basename(infile))[0] + '_' + storm_name
    dataset_name = os.path.splitext(os.path.basename(infile))[0]
    start_time = num2date(time.min(), time_units)
    stop_time = num2date(time.max(), time_units)

    # Show lon/lat
    # import matplotlib.pyplot as plt
    # #plt.plot(lon.flatten(), lat.flatten(), '+b')
    # valid = np.ma.getmaskarray(wind) == False
    # plt.plot(lon[valid], lat[valid], '+b')
    # plt.plot(gcplon.flatten(), gcplat.flatten(), 'ok')
    # gcplon = gcplon.reshape(ngcps)
    # gcplat = gcplat.reshape(ngcps)
    # plt.plot(gcplon[0, :], gcplat[0, :], 'om')
    # plt.plot(gcplon[-1, :], gcplat[-1, :], 'om')
    # plt.plot(gcplon[:, 0], gcplat[:, 0], 'or')
    # plt.plot(gcplon[:, -1], gcplat[:, -1], 'or')
    # xlim = plt.xlim()
    # plt.plot(xlim, [85.05] * 2, '-r')
    # plt.plot(xlim, [-85.05] * 2, '-r')
    # plt.show()
    # import pdb ; pdb.set_trace()

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    metadata['product_name'] = 'SMOSSTORM_AMSR2_wind'
    metadata['name'] = dataset_name
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'JAXA'
    metadata['processing_center'] = 'SOLAB'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'wind speed'
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                gcppix, gcplin)
    band = []
    indndv = np.where(np.ma.getmaskarray(wind))
    offset, scale = vmin, (vmax-vmin)/254.
    _wind = np.clip(np.ma.getdata(wind), vmin, vmax)
    array = np.round((_wind - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('matplotlib_jet',
                                         vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'wind speed', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
