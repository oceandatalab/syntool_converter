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
import pyproj
import gdal
import osr
from pyresample import geometry, kd_tree
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import multiprocessing as mp
import random
import string


def get_gcps_from_bowtie(lon, lat, scansize, ngcps=(10, 10)):
    """
    """
    ysize, xsize = lon.shape
    nscans = ysize / int(scansize)
    indscans = np.round(np.linspace(0, nscans - 1, num=ngcps[0])).astype('int')
    indlines = indscans * int(scansize) + int(scansize) / 2

    gcplon = np.zeros(ngcps)
    gcplat = np.zeros(ngcps)
    gcpnpixel = np.zeros(ngcps)
    gcpnline = (indlines + 0.5 - 0.5) / ysize
    gcpnline = np.tile(gcpnline.reshape(-1, 1), [1, ngcps[1]])

    geod = pyproj.Geod(ellps='WGS84')
    ndistq = np.linspace(0, 1, num=ngcps[1])
    for cnt, indline in enumerate(indlines):
        inv = geod.inv(lon[indline, :], lat[indline, :], lon[indline - 1, :],
                       lat[indline - 1, :])
        forw, dummy, dist = inv
        fwd = geod.fwd(lon[indline, :], lat[indline, :], forw, dist / 2.)
        lonscan, latscan, dummy = fwd
        spacings = geod.inv(lonscan[:-1], latscan[:-1], lonscan[1:],
                            latscan[1:])[2]
        spacing0 = spacings[0] * 3 / 4 - spacings[1] / 4
        spacing1 = spacings[-1] * 3 / 4 - spacings[-2] / 4
        tdist = spacing0 + np.sum(spacings) + spacing1
        ndist = np.hstack((spacing0, spacings)).cumsum() / tdist
        ind2 = [abs(ndist - i).argmin() for i in ndistq]
        gcplon[cnt, :] = lonscan[ind2]
        gcplat[cnt, :] = latscan[ind2]
        gcpnpixel[cnt, :] = ndist[ind2]

    return (gcplon, gcplat, gcpnpixel, gcpnline)


def get_gcps_from_bowtie_old(lon, lat, scansize, ngcps=(10, 10)):
    """
    """
    ysize, xsize = lon.shape
    nscans = ysize / int(scansize)
    midscanlines = np.arange(nscans) * int(scansize) + int(scansize) / 2
    midscanlon = lon[midscanlines.reshape(-1, 1), [[0, xsize / 2, -1]]]
    midscanlat = lat[midscanlines.reshape(-1, 1), [[0, xsize / 2, -1]]]

    if isinstance(midscanlon, ma.MaskedArray) or \
       isinstance(midscanlat, ma.MaskedArray):
        mask = ma.getmaskarray(midscanlon) | ma.getmaskarray(midscanlat)
        valid = np.where(mask.sum(axis=1) == 0)[0]
        if valid.size == nscans:
            pass
        elif valid.size < ngcps[0]:
            raise Exception('Too much invalid lon/lat for setting GCPs !')
        else:
            midscanlon = midscanlon[valid, :]
            midscanlat = midscanlat[valid, :]
            midscanlines = midscanlines[valid]
            nscans = valid.size

    indscan = np.round(np.linspace(0, nscans - 1, num=ngcps[0])).astype('int')
    gcplon = np.zeros(ngcps)
    gcplon[:, 0] = midscanlon[indscan, 0]
    gcplon[:, -1] = midscanlon[indscan, -1]
    gcplat = np.zeros(ngcps)
    gcplat[:, 0] = midscanlat[indscan, 0]
    gcplat[:, -1] = midscanlat[indscan, -1]
    gcpnline = (midscanlines[indscan] + 0.5) / ysize
    gcpnline = np.tile(gcpnline.reshape(-1, 1), [1, ngcps[1]])

    gcpnpixel = np.linspace(0.5 / xsize, 1 - 0.5 / xsize, num=ngcps[1])
    gcpnpixel = np.tile(gcpnpixel.reshape(1, -1), [ngcps[0], 1])
    if ngcps[1] > 2:
        geod = pyproj.Geod(ellps='WGS84')
        for ind in range(ngcps[0]):
            pts = np.array(geod.npts(gcplon[ind, 0], gcplat[ind, 0],
                                     gcplon[ind, -1], gcplat[ind, -1],
                                     ngcps[1] - 2))
            gcplon[ind, 1:-1] = pts[:, 0]
            gcplat[ind, 1:-1] = pts[:, 1]

    return (gcplon, gcplat, gcpnpixel, gcpnline)


def id_generator(size=10, chars=string.ascii_letters + string.digits):
    """
    """
    return ''.join(random.choice(chars) for _ in range(size))


def get_gdal_transformer(gcplon, gcplat, gcppixel, gcpline, xsize, ysize):
    """
    """
    gdal.SetCacheMax(2**30)
    gcps = []
    for i, j, k, l in zip(gcplon.flat, gcplat.flat, gcppixel.flat, gcpline.flat):
        gcps.append(gdal.GCP(float(i), float(j), 0., float(k), float(l)))
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    proj = srs.ExportToWkt()
    drv = gdal.GetDriverByName('MEM')
    dset = drv.Create('tmp', int(xsize), int(ysize))
    #dset = drv.Create(id_generator(), int(xsize), int(ysize))
    dset.SetGCPs(gcps, proj)
    #options = ['']
    #options = ['MAX_GCP_ORDER=3']
    options = ['MAX_GCP_ORDER=-1']
    transformer = gdal.Transformer(dset, None, options)
    return transformer


def transform_points(gcplon, gcplat, gcppixel, gcpline, xsize, ysize,
                     dsttosrc, inputs):
    """
    """
    transformer = get_gdal_transformer(gcplon, gcplat, gcppixel, gcpline,
                                       xsize, ysize)
    outputs = np.array(transformer.TransformPoints(dsttosrc, inputs)[0],
                       dtype='float32')[:, 0:2]
    return outputs


def get_points_from_gcps(gcplon, gcplat, gcppixel, gcpline, xsize, ysize,
                         dsttosrc, lonorpix, latorlin, nprocs=1):
    """
    """
    inputs = np.vstack((lonorpix.flatten(), latorlin.flatten())).transpose()
    ncpus = min(nprocs, mp.cpu_count())
    if ncpus > 1:
        try:
            limits = np.round(np.linspace(0, inputs.shape[0],
                                          num=ncpus + 1)).astype('int32')
            slices = [slice(limits[i], limits[i+1]) for i in range(ncpus)]
            pool = mp.Pool(processes=ncpus)
            jobs = []
            for sli in slices:
                args = (gcplon, gcplat, gcppixel, gcpline, xsize, ysize,
                        dsttosrc, inputs[sli, :])
                jobs.append(pool.apply_async(transform_points, args))
            outputs = np.vstack([job.get() for job in jobs])
            pool.close()
            pool.join()
        except:
            print 'get_points_from_gcps multiprocessing failed, try without it'
            outputs = transform_points(gcplon, gcplat, gcppixel, gcpline,
                                       xsize, ysize, dsttosrc, inputs)
    else:
        outputs = transform_points(gcplon, gcplat, gcppixel, gcpline,
                                   xsize, ysize, dsttosrc, inputs)
    return outputs.transpose().reshape((2,) + lonorpix.shape)


def get_points_from_transformer(transformer, dsttosrc, lonorpix, latorlin):
    """
    """
    inputs = np.vstack((lonorpix.flatten(), latorlin.flatten())).transpose()
    outputs = np.array(transformer.TransformPoints(dsttosrc, inputs)[0],
                       dtype='float32')[:, 0:2]
    return outputs.transpose().reshape((2,) + lonorpix.shape)


# TO DO for optimisation :
# - cut scans into valid areas for triangulation/interpolation
# - multiprocessing
def resample_bowtie_linear(pixels, lines, values, scansize, rspxsize, rspysize,
                           show=False):
    """
    """
    # import time
    # time_start = time.clock()
    # time_triangle = 0.
    # time_interpfunc = 0.
    # time_interp = 0.

    ysize, xsize = pixels.shape
    nscans = ysize / scansize
    valndim = values.ndim
    if valndim == 2:
        values = values.reshape(values.shape + (1,))
        nvar = 1
    else:
        nvar = values.shape[2]
    rspval = np.zeros((rspysize, rspxsize, nvar), dtype='float32')
    rspcnt = np.zeros((rspysize, rspxsize, nvar), dtype='int32')

    for iscan in range(nscans):

        #print iscan
        yind0 = iscan * scansize
        yind1 = yind0 + scansize
        ysl = slice(yind0, yind1)
        xsl = slice(0, xsize)

        # Extend left/right of scan (extrapolation)
        pixlr = 2 * pixels[ysl, [0, -1]] - pixels[ysl, [1, -2]]
        linlr = 2 * lines[ysl, [0, -1]] - lines[ysl, [1, -2]]
        vallr = 2 * values[ysl, [0, -1], :] - values[ysl, [1, -2], :]
        ptsx = np.hstack((pixlr[:, [0]], pixels[ysl, xsl], pixlr[:, [1]]))
        ptsy = np.hstack((linlr[:, [0]], lines[ysl, xsl], linlr[:, [1]]))
        ptsv = np.hstack((vallr[:, [0], :], values[ysl, xsl, :], vallr[:, [1], :]))

        # Extend bottom/up of scan
        # bottom first scan and up last scan : extrapolation
        # others : use last/first line of previous/next scan if outside scan
        #          use NaN otherwise
        # pixbu = 2 * ptsx[[0, -1], :] - ptsx[[1, -2], :]
        # linbu = 2 * ptsy[[0, -1], :] - ptsy[[1, -2], :]
        # valbu = np.empty((2, xsize + 2, nvar))
        # valbu.fill(np.nan)
        # if iscan == 0:
        #     valbu[0, :, :] = 2 * ptsv[0, :, :] - ptsv[1, :, :]
        # else:
        #     inside = np.where(lines[yind0, :] > lines[yind0 - 1, :])[0]
        #     pixbu[0, inside + 1] = pixels[yind0 - 1, inside]
        #     linbu[0, inside + 1] = lines[yind0 - 1, inside]
        #     valbu[0, inside + 1, :] = values[yind0 - 1, inside, :]
        # if iscan == nscans - 1:
        #     valbu[1, :, :] = 2 * ptsv[-1, :, :] - ptsv[-2, :, :]
        # else:
        #     inside = np.where(lines[yind1 - 1, :] < lines[yind1, :])[0]
        #     pixbu[1, inside + 1] = pixels[yind1, inside]
        #     linbu[1, inside + 1] = lines[yind1, inside]
        #     valbu[1, inside + 1, :] = values[yind1, inside, :]
        # ptsx = np.vstack((pixbu[[0], :], ptsx, pixbu[[1], :]))
        # ptsy = np.vstack((linbu[[0], :], ptsy, linbu[[1], :]))
        # ptsv = np.vstack((valbu[[0], :, :], ptsv, valbu[[1], :, :]))

        # Extend bottom/up of scan (extrapolation)
        # And add NaN lines before/after to avoid interpolation far away
        pixbu = 2 * ptsx[[0, -1], :] - ptsx[[1, -2], :]
        linbu = 2 * ptsy[[0, -1], :] - ptsy[[1, -2], :]
        valbu = 2 * ptsv[[0, -1], :, :] - ptsv[[1, -2], :, :]
        pixbunan = 3 * ptsx[[0, -1], :] - 2 * ptsx[[1, -2], :]
        linbunan = 3 * ptsy[[0, -1], :] - 2 * ptsy[[1, -2], :]
        valbunan = np.empty((2, xsize + 2, nvar))
        valbunan.fill(np.nan)
        ptsx = np.vstack((pixbunan[[0], :], pixbu[[0], :], ptsx,
                          pixbu[[1], :], pixbunan[[1], :]))
        ptsy = np.vstack((linbunan[[0], :], linbu[[0], :], ptsy,
                          linbu[[1], :], linbunan[[1], :]))
        ptsv = np.vstack((valbunan[[0], :, :], valbu[[0], :, :], ptsv,
                          valbu[[1], :, :], valbunan[[1], :, :]))

        # Compute min/max output coordinates for interpolation
        indfin = np.where(np.isfinite(ptsv).sum(axis=2) > 0)
        if indfin[0].size == 0:
            continue
        rspx0 = np.clip(np.ceil(ptsx[indfin].min()),
                        0, rspxsize - 1).astype('int32')
        rspx1 = np.clip(np.floor(ptsx[indfin].max()),
                        0, rspxsize - 1).astype('int32')
        if rspx0 > rspx1:
            continue
        rspy0 = np.clip(np.ceil(ptsy[indfin].min()),
                        0, rspysize - 1).astype('int32')
        rspy1 = np.clip(np.floor(ptsy[indfin].max()),
                        0, rspysize - 1).astype('int32')
        if rspy0 > rspy1:
            continue

        # Interpolate
        # time0 = time.clock()
        triangles = Delaunay(np.vstack((ptsx.flatten(),
                                        ptsy.flatten())).transpose())
        # time_triangle += time.clock() - time0
        for ivar in range(nvar):
            # time0 = time.clock()
            interpfunc = LinearNDInterpolator(triangles,
                                              ptsv[:, :, ivar].flatten(),
                                              fill_value=np.nan)
            # time_interpfunc += time.clock() - time0
            # time0 = time.clock()
            tmpval = interpfunc(np.arange(rspx0, rspx1 + 1).reshape((1, -1)),
                                np.arange(rspy0, rspy1 + 1).reshape((-1, 1)))
            # time_interp += time.clock() - time0
            valid = np.where(np.isfinite(tmpval))
            rspval[valid[0] + rspy0, valid[1] + rspx0, ivar] += tmpval[valid]
            rspcnt[valid[0] + rspy0, valid[1] + rspx0, ivar] += 1

    # print 'time_total : {}'.format(time.clock() - time_start)
    # print 'time_triangle : {}'.format(time_triangle)
    # print 'time_interpfunc : {}'.format(time_interpfunc)
    # print 'time_interp : {}'.format(time_interp)

    rspmask = (rspcnt == 0)
    rspval[~rspmask] /= rspcnt[~rspmask]
    rspvalmask = ma.MaskedArray(rspval, rspmask)

    if show == True:
        import matplotlib.pyplot as plt
        for ivar in range(nvar):
            for img in (rspvalmask[:, :, ivar], rspcnt[:, :, ivar]):
                plt.figure()
                plt.imshow(img, interpolation='nearest')
                plt.colorbar()
        plt.show()

    if valndim == 2:
        return rspvalmask[:, :, 0]
    else:
        return rspvalmask


def resample_gauss(lon, lat, var, rsplon, rsplat, radius, sigma,
                   fill_value=None, nprocs=1, show=False):
    """
    """
    swath = geometry.SwathDefinition(lons=lon, lats=lat)
    rspswath = geometry.SwathDefinition(lons=rsplon, lats=rsplat)
    rspvar = kd_tree.resample_gauss(swath, var, rspswath, radius, sigma,
                                    fill_value=fill_value, nprocs=nprocs,
                                    with_uncert=show)
    if show == True:
        import matplotlib.pyplot as plt
        imgs = rspvar
        if len(var.shape) == 2:
            imgs = (var,) + imgs
        for img in imgs:
            plt.figure()
            plt.imshow(img, interpolation='nearest')
            plt.colorbar()
        plt.show()
        rspvar = rspvar[0]
    return rspvar.astype(var.dtype)


def resample_nearest(lon, lat, var, rsplon, rsplat, radius,
                     fill_value=None, nprocs=1, show=False):
    """
    """
    swath = geometry.SwathDefinition(lons=lon, lats=lat)
    rspswath = geometry.SwathDefinition(lons=rsplon, lats=rsplat)
    rspvar = kd_tree.resample_nearest(swath, var, rspswath, radius,
                                      fill_value=fill_value, nprocs=nprocs)
    if show == True:
        import matplotlib.pyplot as plt
        imgs = (rspvar,)
        if len(var.shape) == 2:
            imgs = (var,) + imgs
        for img in imgs:
            plt.figure()
            plt.imshow(img, interpolation='nearest')
            plt.colorbar()
        plt.show()
    return rspvar.astype(var.dtype)
