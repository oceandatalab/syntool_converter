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

# IMPORTANT : use agg backend for fig2rgb() working
# the next two lines must be called before any matplotlib import
import matplotlib
matplotlib.use('agg')
#matplotlib.use('ps') # does not support tostring_rgb()
# matplotlib.rcParams['lines.antialiased'] = False
# matplotlib.rcParams['patch.antialiased'] = False
# matplotlib.rcParams['text.antialiased'] = False
from sar.utils.factory import sarimage
from sar.transform.sarimage2sarxspec import sarimage2sarxspec_loop
from sar.processing.make_sarxspec import make_sarxspec_fig
from sar.render.palette import getColorMap
import syntool_converter.utils.syntoolformat as stfmt
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import pyproj
import re


def fig2rgb(fig):
    """
    """
    fig.canvas.draw()
    #fig.canvas.print_png('/tmp/test_sm_xspec.png')
    wid, hei = fig.canvas.get_width_height()
    rgb = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
    rgb.shape = (hei, wid, 3)
    return rgb


def nuniqcolors(rgb):
    """
    """
    col = rgb[:, :, 0]*65536 + rgb[:, :, 1]*256 + rgb[:, :, 2]
    uniqcol = np.unique(col)
    return uniqcol.size


def fig2imgpal(fig, cmap):
    """
    """
    # Convert fig to PIL
    from PIL import Image
    fig.canvas.draw()
    wid, hei = fig.canvas.get_width_height()
    img = Image.fromstring('RGB', [wid, hei], fig.canvas.tostring_rgb())
    # Make PIL palette from cmap
    ncolors = np.minimum(cmap.N, 254)
    colors = cmap(np.linspace(0, 1, num=ncolors))
    palette = []
    for rgba in colors:
        entry = [int(col*255) for col in rgba[0:-1]]
        palette.extend(entry)
    nfill = 256-ncolors-1
    for i in range(0, nfill):
        palette.extend([int(i*255./(nfill-1.)) for x in range(3)])
    palette.extend([0, 0, 0]) # for nodatavalues
    # Create a dummy image with PIL palette
    dummy = Image.new('P', [50, 50])
    dummy.putpalette(palette)
    # Quantize
    palimg = img.quantize(palette=dummy)
    if palimg.getextrema()[1] == 255:
        raise Exception('Unexpected quantization.')
    #print palimg.getextrema()
    return (np.asarray(palimg), palimg.getpalette())


def palette2colortable(pal):
    """
    """
    import gdal
    colortable = gdal.ColorTable()
    for index in range(0, len(pal)/3):
        entry = [int(col) for col in pal[index*3:index*3+3]]
        colortable.SetColorEntry(index, tuple(entry))
    return colortable


def get_cmaps(ncolors=None):
    """
    """
    cmap_re = getColorMap(rgbFile='wind.pal')
    cmap_im = plt.cm.PuOr
    cmaps = [cmap_re, cmap_im]
    if ncolors != None:
        for i in [0, 1]:
            colors = cmaps[i](np.linspace(0, 1, num=min([cmaps[i].N, ncolors])))
            cmaps[i] = LinearSegmentedColormap.from_list('custom', colors, N=ncolors)
    return cmaps


def sar_xspec(infile, outdir, pngkml=False, vmax_re=None, vmax_im=None,
              make_rgb=True, ncolors=74):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    sarim = sarimage(infile)
    mission = sarim.get_info('mission')
    if mission == 'S1A':
        sensor_name = 'Sentinel-1A'
        sensor_platform = 'Sentinel-1A'
        source_provider = 'ESA'
    else:
        raise Exception('S1A mission expected.')
    product = sarim.get_info('product')
    if product != 'SLC':
        raise Exception('SLC expected.')
    timefmt = '%Y-%m-%dT%H:%M:%S.%f'
    start_time = datetime.strptime(sarim.get_info('start_time'), timefmt)
    stop_time = datetime.strptime(sarim.get_info('stop_time'), timefmt)
    sensor_pass = sarim.get_info('pass')
    sensor_mode = sarim.get_info('mode')
    sensor_swath = sarim.get_info('swath')
    sensor_polarisation = sarim.get_info('polarisation')
    datagroup = sarim.get_info('safe_name').replace('.SAFE', '')
    pid = datagroup.split('_')[-1]
    dataname = os.path.splitext(os.path.basename(infile))[0] + '-' + pid
    # Compute SAR Xspec and make figures
    if sensor_mode == 'WV':
        azi_periodo_size = 1024
        azi_dist, ran_dist = 20000., 20000. # ignored in WV case
        xspec_size = (512, 512)
        fontsize = 'small'
    elif re.match(r'^S[1-6]$', sensor_mode) != None:
        azi_periodo_size = 1024
        azi_dist, ran_dist = 10000., 10000.
        xspec_size = (512, 512) #(256, 256)
        fontsize = 'small' #'x-small'
    elif sensor_mode in ['IW', 'EW']:
        azi_periodo_size = 512
        azi_dist, ran_dist = 10000., 10000.
        xspec_size = (512, 512) #(256, 256)
        fontsize = 'small' #'x-small'
    else:
        raise Exception('Which settings for this mode ?')
    sarxspec = sarimage2sarxspec_loop(sarim, azi_dist=azi_dist, ran_dist=ran_dist,
                                      azi_periodo_size=azi_periodo_size)
    cmap_re, cmap_im = get_cmaps(ncolors=ncolors)
    fig_re = make_sarxspec_fig(sarxspec, part='real', tau=1, kmax=2*np.pi/75,
                               kmin=2*np.pi/400, xspec_size=xspec_size,
                               uniq_vmax=True, north_oriented=True,
                               klim=[2*np.pi/400, 2*np.pi/200, 2*np.pi/100],
                               north_arrow=False, index_pos=None, vmax_pos='tr',
                               nvar_pos=None, fontsize=fontsize,
                               vmax=vmax_re, cmap=cmap_re)
    if sensor_mode == 'WV':
        ax = fig_re.gca()
        imnum = sarim.get_info('image_number')
        imnumtxt = '#{:03d}'.format(imnum)
        ax.text(0.51, 0.99, imnumtxt, transform=ax.transAxes,
                ha='left', va='top', fontsize=fontsize)
    if make_rgb == True:
        rgb_re = fig2rgb(fig_re)
        #print nuniqcolors(rgb_re)
    else:
        img_re, pal_re = fig2imgpal(fig_re, cmap_re)
    plt.close(fig_re)
    fig_im = make_sarxspec_fig(sarxspec, part='imag', tau=1, kmax=2*np.pi/75,
                               kmin=2*np.pi/400, xspec_size=xspec_size,
                               uniq_vmax=True, north_oriented=True,
                               klim=[2*np.pi/400, 2*np.pi/200, 2*np.pi/100],
                               north_arrow=False, index_pos=None, vmax_pos='tr',
                               nvar_pos=None, fontsize=fontsize,
                               vmax=vmax_im, cmap=cmap_im)
    if sensor_mode == 'WV':
        ax = fig_im.gca()
        imnum = sarim.get_info('image_number')
        imnumtxt = '#{:03d}'.format(imnum)
        ax.text(0.51, 0.99, imnumtxt, transform=ax.transAxes,
                ha='left', va='top', fontsize=fontsize)
    if make_rgb == True:
        rgb_im = fig2rgb(fig_im)
        #print nuniqcolors(rgb_im)
    else:
        img_im, pal_im = fig2imgpal(fig_im, cmap_im)
    plt.close(fig_im)
    if make_rgb == True:
        nlin, npix = rgb_re.shape[0:2]
    else:
        nlin, npix = img_re.shape
    # Handle GCPS
    # geoloc = sarim.get_info('geolocation_grid')
    # pix = np.array([0, geoloc['npixels']-1, geoloc['npixels']-1, 0])
    # lin = np.array([0, 0, geoloc['nlines']-1, geoloc['nlines']-1])
    # gcplon = geoloc['longitude'][lin, pix]
    # gcplat = geoloc['latitude'][lin, pix]
    # gcphei = np.zeros(4)
    # gcppix = np.array([0, 512, 512, 0])
    # gcplin = np.array([0, 0, 512, 512])
    #############################################
    # geoloc = sarim.get_info('geolocation_grid')
    # gcplon = geoloc['longitude']
    # gcplat = geoloc['latitude']
    # gcphei = np.zeros(gcplon.shape)
    # geod = pyproj.Geod(ellps='WGS84')
    # nglin, ngpix = geoloc['nlines'], geoloc['npixels']
    # ra_geo_spacing = geod.inv(gcplon[nglin/2, 0:-1], gcplat[nglin/2, 0:-1],
    #                           gcplon[nglin/2, 1:], gcplat[nglin/2, 1:])[2]
    # ra_geo_dist = np.hstack((0., ra_geo_spacing.cumsum()))
    # ra_geo_ndist = ra_geo_dist/ra_geo_dist[-1]
    # gcppix = np.tile((ra_geo_ndist*npix).reshape((1, -1)), (nglin, 1))
    # az_geo_spacing = geod.inv(gcplon[0:-1, ngpix/2], gcplat[0:-1, ngpix/2],
    #                           gcplon[1:, ngpix/2], gcplat[1:, ngpix/2])[2]
    # az_geo_dist = np.hstack((0., az_geo_spacing.cumsum()))
    # az_geo_ndist = az_geo_dist/az_geo_dist[-1]
    # gcplin = np.tile((az_geo_ndist*nlin).reshape((-1, 1)), (1, ngpix))
    # import pdb ; pdb.set_trace()
    #############################################
    #import pdb ; pdb.set_trace()
    ext_min = sarxspec[0][0].get_info('extent')[0:2]
    ext_max = sarxspec[-1][-1].get_info('extent')[2:4]
    # geoloc = sarim.get_info('geolocation_grid')
    # nglin, ngpix = geoloc['nlines'], geoloc['npixels']
    nglin, ngpix = len(sarxspec)+1, len(sarxspec[0])+1
    pix = np.hstack((np.round(np.linspace(ext_min[1], ext_max[1], num=ngpix)),
                     np.ones(nglin)*ext_max[1],
                     np.round(np.linspace(ext_max[1], ext_min[1], num=ngpix)),
                     np.ones(nglin)*ext_min[1]))
    lin = np.hstack((np.ones(ngpix)*ext_min[0],
                     np.round(np.linspace(ext_min[0], ext_max[0], num=nglin)),
                     np.ones(ngpix)*ext_max[0],
                     np.round(np.linspace(ext_max[0], ext_min[0], num=nglin))))
    lon, lat = np.zeros(pix.size), np.zeros(pix.size)
    for ipt in range(pix.size):
        ext = [lin[ipt], pix[ipt], lin[ipt], pix[ipt]]
        lon[ipt] = sarim.get_data('lon', extent=ext, spacing=1)
        lat[ipt] = sarim.get_data('lat', extent=ext, spacing=1)
    ndist = np.zeros(pix.size)
    lim = [0, ngpix, ngpix+nglin, 2*ngpix+nglin, 2*ngpix+2*nglin]
    geod = pyproj.Geod(ellps='WGS84')
    for iside in range(4):
        pt0, pt1 = lim[iside], lim[iside+1]-1
        ddist = geod.inv(lon[pt0:pt1], lat[pt0:pt1], lon[pt0+1:pt1+1],
                         lat[pt0+1:pt1+1])[2]
        dist = ddist.cumsum()
        ndist[pt0:pt1+1] = np.hstack((0., dist))/dist.max()
    gcppix = np.hstack((ndist[lim[0]:lim[1]-1]*npix, np.ones(nglin-1)*npix,
                        (1-ndist[lim[2]:lim[3]-1])*npix, np.zeros(nglin-1)))
    gcplin = np.hstack((np.zeros(ngpix-1), ndist[lim[1]:lim[2]-1]*nlin,
                        np.ones(ngpix-1)*nlin, (1-ndist[lim[3]:lim[4]-1])*nlin))
    gcplon = np.hstack((lon[lim[0]:lim[1]-1], lon[lim[1]:lim[2]-1],
                        lon[lim[2]:lim[3]-1], lon[lim[3]:lim[4]-1]))
    gcplat = np.hstack((lat[lim[0]:lim[1]-1], lat[lim[1]:lim[2]-1],
                        lat[lim[2]:lim[3]-1], lat[lim[3]:lim[4]-1]))
    gcphei = np.zeros(gcplon.size)
    #import pdb ; pdb.set_trace()
    if gcplon.min() < -135 and gcplon.max() > 135:
        gcplon[np.where(gcplon < 0)] += 360.
    #############################################
    if sensor_pass == 'Descending':
        gcppix = npix-gcppix
        gcplin = nlin-gcplin
    gcplin = nlin-gcplin # because fig will be read and wrote from top to bottom
    # Loop on part and write
    for part in ['real', 'imag']:
        print part
        if part == 'real':
            product = 'SAR_cross-spectrum_real'
            nameext = '-xspec_re'
            if make_rgb == True:
                rgb = rgb_re
            else:
                img = img_re
                pal = pal_re
        elif part == 'imag':
            product = 'SAR_cross-spectrum_imaginary'
            nameext = '-xspec_im'
            if make_rgb == True:
                rgb = rgb_im
            else:
                img = img_im
                pal = pal_im
        # Construct metadata/geolocation/band(s)
        print 'Construct metadata/geolocation/band(s)'
        metadata = {}
        (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                          units='ms')
        metadata['product_name'] = product
        metadata['name'] = dataname + nameext
        metadata['datetime'] = dtime
        metadata['time_range'] = time_range
        metadata['source_URI'] = infile
        metadata['source_provider'] = source_provider
        metadata['processing_center'] = 'OceanDataLab'
        metadata['conversion_software'] = 'Syntool'
        metadata['conversion_version'] = '0.0.0'
        metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
        metadata['parameter'] = ''
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
        if make_rgb == True:
            band.append({'array': rgb[:, :, 0]})
            band.append({'array': rgb[:, :, 1]})
            band.append({'array': rgb[:, :, 2]})
        else:
            band.append({'array': img, 'nodatavalue': 255,
                         'colortable': palette2colortable(pal)})
        # Write geotiff
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
        # Write projected png/kml
        if pngkml == True:
            print 'Write projected png/kml'
            stfmt.write_pngkml_proj(tifffile)
