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
from cerbere.mapper.safeocnncfile import SAFEOCNNCFile
import os
import numpy as np
import matplotlib.pyplot as plt
from sar.render.palette import getColorMap
import pyproj
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
from scipy.stats import scoreatpercentile


def ax_txt_corner(ax, txt, pos, **kwargs):
    """
    """
    pos2set = {'tl': {'x':0.01, 'y':0.99, 'ha':'left', 'va':'top'},
               'tr': {'x':0.99, 'y':0.99, 'ha':'right', 'va':'top'},
               'bl': {'x':0.01, 'y':0.01, 'ha':'left', 'va':'bottom'},
               'br': {'x':0.99, 'y':0.01, 'ha':'right', 'va':'bottom'},
               'tcr': {'x':0.51, 'y':0.99, 'ha':'left', 'va':'top'}}
    xlim = ax.get_xlim()
    x = pos2set[pos]['x']
    xpos = xlim[0] + x*(xlim[1]-xlim[0])
    ylim = ax.get_ylim()
    y = pos2set[pos]['y']
    ypos = ylim[0] + y*(ylim[1]-ylim[0])
    ha = pos2set[pos]['ha']
    va = pos2set[pos]['va']
    ax.text(xpos, ypos, txt, ha=ha, va=va, **kwargs)


def make_spec_fig(spec, kctr, phictr, heading,
                  npartitions, partitions, hs, flag, imnum=None,
                  spec_size=(512, 512), kmax=2*np.pi/75,
                  kcircle=[2*np.pi/400, 2*np.pi/200, 2*np.pi/100],
                  pcolor=['r', 'b', 'g', 'c'], ppos=['tl', 'tr', 'br', 'bl'],
                  cmap=None, fontsize=None):
    """
    """
    # Constants
    naz, nra, nphi, nk = spec.shape
    klim = np.concatenate([[kctr[0] - (kctr[1] - kctr[0]) / 2.],
                           (kctr[0:-1] + kctr[1:]) / 2.,
                           [kctr[-1] + (kctr[-1] - kctr[-2]) / 2.]])
    dphi = 360./phictr.size
    philim = np.append(phictr - dphi/2, phictr[-1] + dphi/2)
    klimgrid, philimgrid = np.meshgrid(klim, philim)
    theta = 2*np.pi*np.linspace(0, 1, num=361)
    if nra == 1 and naz == 1:
        vmax = spec.max()
    else:
        indsea = np.where(flag == 0)
        if indsea[0].size == 0:
            vmax = spec.max()
        else:
            #vmax = scoreatpercentile(spec[indsea[0], indsea[1], :, :], 99.9)
            seaspec = spec[indsea[0], indsea[1], :, :]
            vmaxs = seaspec.reshape((seaspec.shape[0], -1)).max(axis=1)
            vmax = vmaxs.mean()

    # Create figure
    dpi = 100.
    fig_size = np.array((nra, naz), dtype='float')*spec_size
    fig_size_inch = fig_size/dpi
    fig = plt.figure(figsize=fig_size_inch, dpi=dpi, facecolor='w')

    # Make figure
    for iaz in range(naz):
        for ira in range(nra):

            # Set xlim, ylim for pcolormesh
            #philimgrid2 = philimgrid
            philimgrid2 = 90. - philimgrid + heading[iaz, ira]
            if np.sin((90 - heading[iaz, ira]) * np.pi / 180) < 0:
                # descending pass
                philimgrid2 += 180
            xlim = klimgrid * np.cos(philimgrid2 * np.pi / 180)
            ylim = klimgrid * np.sin(philimgrid2 * np.pi / 180)

            # Create axes and plot
            if np.sin((90 - heading[iaz, ira]) * np.pi / 180) < 0:
                # descending pass
                left = (nra-1-ira)*spec_size[0]/fig_size[0]
                bottom = (naz-1-iaz)*spec_size[1]/fig_size[1]
            else:
                left = ira*spec_size[0]/fig_size[0]
                bottom = iaz*spec_size[1]/fig_size[1]
            rect = [left, bottom, spec_size[0]/fig_size[0],
                    spec_size[1]/fig_size[1]]
            ax = fig.add_axes(rect)
            axkmin = -kmax
            axkmax = kmax
            ax.set_ylim([axkmin, axkmax])
            ax.set_xlim([axkmin, axkmax])

            ax.pcolormesh(xlim, ylim, spec[iaz, ira, :, :],
                          cmap=cmap, vmin=0., vmax=vmax)

            for partnum in range(0, npartitions):
                indpart = np.where(partitions[iaz, ira, :, :] == partnum)
                for x, y in zip(indpart[0], indpart[1]):
                    if y == 0:
                        plt.plot(xlim[x:x+2, y], ylim[x:x+2, y],
                                 color=pcolor[partnum])
                    else:
                        partnumneigh = partitions[iaz, ira, x, y-1]
                        if partnumneigh >= npartitions:
                            plt.plot(xlim[x:x+2, y], ylim[x:x+2, y],
                                     color=pcolor[partnum])
                        elif partnumneigh != partnum:
                            plt.plot(xlim[x:x+2, y], ylim[x:x+2, y],
                                     color=pcolor[partnum], alpha=0.5)
                    if y == nk-1:
                        plt.plot(xlim[x:x+2, y+1], ylim[x:x+2, y+1],
                                 color=pcolor[partnum])
                    else:
                        partnumneigh = partitions[iaz, ira, x, y+1]
                        if partnumneigh >= npartitions:
                            plt.plot(xlim[x:x+2, y+1], ylim[x:x+2, y+1],
                                     color=pcolor[partnum])
                        elif partnumneigh != partnum:
                            plt.plot(xlim[x:x+2, y+1], ylim[x:x+2, y+1],
                                     color=pcolor[partnum], alpha=0.5)
                    partnumneigh = partitions[iaz, ira, x-1, y]
                    if partnumneigh >= npartitions:
                        plt.plot(xlim[x, y:y+2], ylim[x, y:y+2],
                                 color=pcolor[partnum])
                    elif partnumneigh != partnum:
                        plt.plot(xlim[x, y:y+2], ylim[x, y:y+2],
                                 color=pcolor[partnum], alpha=0.5)
                    partnumneigh = partitions[iaz, ira, (x+1) % nphi, y]
                    if partnumneigh >= npartitions:
                        plt.plot(xlim[x+1, y:y+2], ylim[x+1, y:y+2],
                                 color=pcolor[partnum])
                    elif partnumneigh != partnum:
                        plt.plot(xlim[x+1, y:y+2], ylim[x+1, y:y+2],
                                 color=pcolor[partnum], alpha=0.5)
                hstxt = 'Hs = {0:.2f}m'.format(hs[iaz, ira, partnum])
                ax_txt_corner(ax, hstxt, ppos[partnum],
                              color=pcolor[partnum], fontsize=fontsize)

            if flag[iaz, ira] != 0:
                ax.text(0, 0, 'LAND', ha='center', va='center',
                        color='r', fontsize=fontsize)

            ax.set_axis_off()
            ax.plot([axkmin, axkmax], [0, 0], ':k')
            ax.plot([0, 0], [axkmin, axkmax], ':k')
            if kcircle != None:
                for kcir in kcircle:
                    ax.plot(kcir*np.cos(theta), kcir*np.sin(theta), ':k')
                    kcirstr = '%im' % (np.round(2*np.pi/kcir))
                    ax.text(0, -kcir, kcirstr, ha='center', va='top',
                            fontsize=fontsize)

    if imnum is not None:
        imnumtxt = '#{:03d}'.format(imnum)
        ax_txt_corner(ax, imnumtxt, 'tcr', fontsize=fontsize)
    # plt.show()
    return fig


def fig2rgb(fig):
    """
    """
    fig.canvas.draw()
    #fig.canvas.print_png('/tmp/test_sm_xspec.png')
    wid, hei = fig.canvas.get_width_height()
    rgb = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
    rgb.shape = (hei, wid, 3)
    return rgb


def sar_wave(infile, outdir, pngkml=False):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    sarwave = SAFEOCNNCFile(infile, product='WAVE')
    mission = sarwave.read_global_attribute('missionName')
    if mission == 'S1A':
        sensor_name = 'Sentinel-1A'
        sensor_platform = 'Sentinel-1A'
        source_provider = 'ESA'
    else:
        raise Exception('S1A mission expected.')
    # start_time = sarwave.get_start_time() # WARNING : whole SAFE for imagettes !
    # stop_time = sarwave.get_end_time() # WARNING : whole SAFE for imagettes !
    start_t = sarwave.read_global_attribute('firstMeasurementTime')
    if '.' in start_t:
        start_time = datetime.strptime(start_t, '%Y-%m-%dT%H:%M:%S.%fZ')
    else:
        start_time = datetime.strptime(start_t, '%Y-%m-%dT%H:%M:%SZ')
    stop_t = sarwave.read_global_attribute('lastMeasurementTime')
    if '.' in stop_t:
        stop_time = datetime.strptime(stop_t, '%Y-%m-%dT%H:%M:%S.%fZ')
    else:
        stop_time = datetime.strptime(stop_t, '%Y-%m-%dT%H:%M:%SZ')
    heading = sarwave.read_values('oswHeading')
    if np.sin((90 - heading[0, 0]) * np.pi / 180) > 0:
        sensor_pass = 'Ascending'
    else:
        sensor_pass = 'Descending'
    safe_name = os.path.basename(os.path.dirname(os.path.dirname(infile)))
    sensor_mode = safe_name.split('_')[1]
    if sensor_mode not in ['WV', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
        raise Exception('WV/S[1-6] modes expected.')
    sensor_swath = os.path.basename(infile).split('-')[1].upper()
    sensor_polarisation = sarwave.read_global_attribute('polarisation')
    datagroup = safe_name.replace('.SAFE', '')
    pid = datagroup.split('_')[-1]
    dataname = os.path.splitext(os.path.basename(infile))[0] + '-' + pid

    # Make wave spectrum figure
    spec = sarwave.read_values('oswPolSpec')
    k = sarwave.read_values('oswK')
    phi = sarwave.read_values('oswPhi')
    npartitions = sarwave.get_dimsize('oswPartitions')
    partitions = sarwave.read_values('oswPartitions')
    # TMP Bug : there are now 3 partitions, numbered 0, 1 and 3 ...
    if npartitions == 3:
        indp2 = np.where(partitions == 2)
        indp3 = np.where(partitions == 3)
        if indp2[0].size == 0 and indp3[0].size != 0:
            partitions[indp3] = 2
    # /TMP
    hs = sarwave.read_values('oswHs')
    flag = sarwave.read_values('oswLandFlag')
    if sensor_mode == 'WV':
        imnum = int(os.path.splitext(os.path.basename(infile))[0].split('-')[-1])
    else:
        imnum = None
    spec_size = (512, 512)
    fontsize = 'small'
    cmap = getColorMap(rgbFile='wind.pal')
    fig = make_spec_fig(spec, k, phi, heading,
                        npartitions, partitions, hs, flag, imnum=imnum,
                        spec_size=spec_size, fontsize=fontsize, cmap=cmap)
    rgb = fig2rgb(fig)
    plt.close(fig)

    # Make geolocation
    if sensor_mode == 'WV':
        lon = sarwave.read_values('lon')[0, 0]
        lat = sarwave.read_values('lat')[0, 0]
        grdrasize = sarwave.read_values('oswGroundRngSize')[0, 0]
        grdazsize = sarwave.read_values('oswAziSize')[0, 0]
        geod = pyproj.Geod(ellps='WGS84')
        lons = np.repeat(lon, 2)
        lats = np.repeat(lat, 2)
        az = heading[0, 0] + [0., 180.]
        dist = np.repeat(grdazsize / 2., 2)
        lonsmid, latsmid, dummy = geod.fwd(lons, lats, az, dist)
        lons = np.repeat(lonsmid, 2)
        lats = np.repeat(latsmid, 2)
        az = heading[0, 0] + [-90, 90., 90., -90.]
        dist = np.repeat(grdrasize / 2., 4)
        gcplon, gcplat, dummy = geod.fwd(lons, lats, az, dist)
        gcppix = np.array([0, spec_size[0], spec_size[0], 0])
        gcplin = np.array([0, 0, spec_size[1], spec_size[1]])
        if np.sin((90 - heading[0, 0]) * np.pi / 180) < 0:
            # descending pass
            gcppix = spec_size[0] - gcppix
            gcplin = spec_size[1] - gcplin
        gcphei = np.zeros(gcplin.size)
    else:
        gcplon = sarwave.read_values('lon')
        gcplat = sarwave.read_values('lat')
        gcphei = np.zeros(gcplon.shape)
        nra = sarwave.get_dimsize('oswRaSize')
        pix = np.arange(nra) * spec_size[0] + spec_size[0] / 2.
        naz = sarwave.get_dimsize('oswAzSize')
        lin = np.arange(naz-1, -1, -1) * spec_size[1] + spec_size[1] / 2.
        gcppix, gcplin = np.meshgrid(pix, lin)
        if np.sin((90 - heading[0, 0]) * np.pi / 180) < 0:
            # descending pass
            gcppix = nra * spec_size[0] - gcppix
            gcplin = naz * spec_size[1] - gcplin
    if gcplon.min() < -135 and gcplon.max() > 135:
        gcplon[np.where(gcplon < 0)] += 360.

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                          units='ms')
    metadata['product_name'] = 'SAR_wave'
    metadata['name'] = dataname
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = source_provider
    metadata['processing_center'] = '' #'OceanDataLab'
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
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                gcppix, gcplin)
    band = []
    band.append({'array': rgb[:, :, 0]})
    band.append({'array': rgb[:, :, 1]})
    band.append({'array': rgb[:, :, 2]})

    # Write geotiff
    print 'Write geotiff'
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    # Write projected png/kml
    if pngkml == True:
        print 'Write projected png/kml'
        stfmt.write_pngkml_proj(tifffile)
