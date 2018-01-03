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
from netCDF4 import Dataset, num2date
import pyproj
import syntool_converter.utils.syntoolformat as stfmt
import os
import re
import sys
import tools_for_gcp
from scipy import stats, interpolate


def make_gcps_v1(lon, lat, dist_gcp=None):
    """
    """
    if dist_gcp is None:
        dist_gcp = 25.
    geod = pyproj.Geod(ellps='WGS84')
    _, _, dist = geod.inv(lon[:-1], lat[:-1], lon[1:], lat[1:])
    ngcplin = np.ceil(dist.sum() / 1000. / dist_gcp).astype('int') + 1
    _gcpind = np.round(np.linspace(0, lon.size - 1, num=ngcplin)).astype('int')
    _gcppix = np.array([-1., 2.])
    ngcppix = _gcppix.size
    gcppix = np.tile(_gcppix[np.newaxis, :], (ngcplin, 1))
    gcpind = np.tile(_gcpind[:, np.newaxis], (1, ngcppix))
    gcplin = gcpind + 0.5
    swathwid = 5. / 111.32
    _ind = np.minimum(gcpind, lon.size - 2)
    swathdir = np.arctan2(lat[_ind + 1] - lat[_ind], lon[_ind + 1] - lon[_ind])
    gcplon = lon[gcpind] + (gcppix - 0.5)*swathwid*np.cos(swathdir - np.pi/2.)
    gcplat = lat[gcpind] + (gcppix - 0.5)*swathwid*np.sin(swathdir - np.pi/2.)
    gcphei = np.zeros(gcppix.shape)
    return (gcplon, gcplat, gcphei, gcppix, gcplin)


def make_gcps_v2(lon, lat, dist_gcp=None):
    """
    """
    if dist_gcp is None:
        dist_gcp = 1.
    dists = np.sqrt((lon[1:] - lon[:-1])**2. + (lat[1:] - lat[:-1])**2.)
    distscum = np.concatenate(([0.], dists.cumsum()))
    dist = distscum[-1]
    if dist == 0.:
        print('WARNING : null distance in GCPs computation.')
        ngcplin = 2
    else:
        ngcplin = int(np.ceil(dist / dist_gcp)) + 1
    gcpindfunc = interpolate.interp1d(distscum, np.arange(lon.size),
                                      kind='quadratic')
    _gcpind = np.round(gcpindfunc(np.linspace(0., dist,
                       num=ngcplin))).astype('int')
    _gcppix = np.array([-25., 0.5, 26.])
    ngcppix = _gcppix.size
    gcppix = np.tile(_gcppix[np.newaxis, :], (ngcplin, 1))
    gcpind = np.tile(_gcpind[:, np.newaxis], (1, ngcppix))
    gcplin = gcpind + 0.5
    swathwid = 5. / 111.32
    _ind0 = np.minimum(gcpind, lon.size - 2)
    _ind1 = _ind0 + 1
    indsame = np.where((lon[_ind0] == lon[_ind1]) & (lat[_ind0] == lat[_ind1]))
    for igl, igp in zip(indsame[0], indsame[1]):
        if _ind1[igl, igp] < lon.size - 1:
            _ind1[igl, igp] += 1
        else:
            _ind0[igl, igp] -= 1
    if np.min(lon) < -180. :
        lon += 360.
    if np.max(lon) > 360. :
        lon -= 360.
    swathdir = np.arctan2(lat[_ind1] - lat[_ind0], lon[_ind1] - lon[_ind0])
    gcplon = lon[gcpind] + (gcppix - 0.5)*swathwid*np.cos(swathdir - np.pi/2.)
    gcplat = lat[gcpind] + (gcppix - 0.5)*swathwid*np.sin(swathdir - np.pi/2.)
    gcphei = np.zeros(gcppix.shape)
    return (gcplon, gcplat, gcphei, gcppix, gcplin)


def check_gcps(tifffile, lon, lat, gcplon, gcplat):
    """
    """
    from osgeo import gdal
    import matplotlib.pyplot as plt
    import pdb
    print('check GCPs')
    dset = gdal.Open(tifffile)
    tfr = gdal.Transformer(dset, None, ['MAX_GCP_ORDER=-1'])
    #
    pix = np.linspace(0., 1., num=5)
    lin = np.arange(dset.RasterYSize) + 0.5
    xy_grid = np.array((np.tile(pix.reshape((-1, 1)), (1, lin.size)),
                        np.tile(lin.reshape((1, -1)), (pix.size, 1))))
    xy_dims = xy_grid.shape[1:3]
    xy_grid = xy_grid.reshape((2, -1)).transpose()
    _lonlat = np.array(tfr.TransformPoints(0, xy_grid)[0])
    _lon = _lonlat[:, 0].reshape(xy_dims).transpose()
    _lat = _lonlat[:, 1].reshape(xy_dims).transpose()
    plt.figure()
    plt.plot(lon, lat, '-ob', markersize=9)
    plt.plot(_lon.flatten(), _lat.flatten(), 'xg', markersize=9)
    plt.plot(gcplon, gcplat, '+r', markersize=9)
    plt.xlabel('lon')
    plt.ylabel('lat')
    #
    geod = pyproj.Geod(ellps='WGS84')
    _, _, distm = geod.inv(lon, lat, _lon[:, 2], _lat[:, 2])
    distd = np.sqrt((lon - _lon[:, 2])**2. + (lat - _lat[:, 2])**2.)
    wid = np.sqrt((_lon[:, -1] - _lon[:, 0])**2.
                  + (_lat[:, -1] - _lat[:, 0])**2.)
    _wid = 5. / 111.32
    widabserr = np.abs(wid - _wid)
    widrelerr = widabserr / _wid * 100.
    print(gcplon.shape)
    print(distm.min(), distm.max(), distm.mean())
    print(distd.min(), distd.max(), distd.mean())
    print(widabserr.min(), widabserr.max(), widabserr.mean())
    print(widrelerr.min(), widrelerr.max(), widrelerr.mean())
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(lat, distm, '+-')
    plt.xlabel('lat')
    plt.ylabel('dist [m]')
    plt.subplot(2, 2, 2)
    plt.plot(lat, distd, '+-')
    plt.xlabel('lat')
    plt.ylabel('dist [deg]')
    plt.subplot(2, 2, 3)
    plt.plot(lat, wid, '+-')
    plt.plot([lat.min(), lat.max()], [_wid] * 2, '-r')
    plt.xlabel('lat')
    plt.ylabel('width [deg]')
    plt.subplot(2, 2, 4)
    plt.plot(lat, widrelerr, '+-')
    plt.xlabel('lat')
    plt.ylabel('width rel error [%]')
    #
    glon = np.linspace(lon.min() - 1., lon.max() + 1., num=500)
    glat = np.linspace(lat.min() - 1., lat.max() + 1., num=500)
    ll_grid = np.array((np.tile(glon.reshape((-1, 1)), (1, glat.size)),
                        np.tile(glat.reshape((1, -1)), (glon.size, 1))))
    ll_dims = ll_grid.shape[1:3]
    ll_grid = ll_grid.reshape((2, -1)).transpose()
    pixlin = np.array(tfr.TransformPoints(1, ll_grid)[0])
    pix = pixlin[:, 0].reshape(ll_dims).transpose()
    lin = pixlin[:, 1].reshape(ll_dims).transpose()
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.imshow(pix, origin='lower', interpolation='nearest',
               extent=[glon.min(), glon.max(), glat.min(), glat.max()])
    plt.colorbar(label='pixel')
    plt.xlabel('lon')
    plt.ylabel('lat')
    plt.subplot(1, 2, 2)
    plt.imshow(lin, origin='lower', interpolation='nearest',
               extent=[glon.min(), glon.max(), glat.min(), glat.max()])
    plt.colorbar(label='line')
    plt.xlabel('lon')
    plt.ylabel('lat')
    #
    plt.show()
    pdb.set_trace()
