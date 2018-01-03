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

from cerbere.mapper.safemsil1cstitchedfile import safemsil1c_stitching_groups, \
    SAFEMSIL1CStitchedFile
import numpy as np
from scipy.misc import bytescale
import syntool_converter.utils.syntoolformat as stfmt
import os
from datetime import datetime
from osgeo import gdal, ogr, osr
from scipy.ndimage.morphology import binary_dilation


def compute_landmask(mapper, landmaskpath, downsampling=1):
    """
    """
    # Set grid georeferencing
    srs = osr.SpatialReference()
    srs.SetFromUserInput(mapper.read_global_attribute('horizontal_cs_code'))
    geotf = [mapper.read_global_attribute('ulx'),
             mapper.read_global_attribute('xdim') * downsampling,
             0,
             mapper.read_global_attribute('uly'),
             0,
             mapper.read_global_attribute('ydim') * downsampling]
    ysize = int(mapper.read_global_attribute('ncols') / float(downsampling))
    xsize = int(mapper.read_global_attribute('nrows') / float(downsampling))
    # Get landmask layer (polygons) and set spatial filter
    datasource = ogr.Open(landmaskpath)
    layer = datasource.GetLayer(0)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(geotf[0] + geotf[1] * 0, geotf[3] + geotf[5] * 0)
    ring.AddPoint(geotf[0] + geotf[1] * xsize, geotf[3] + geotf[5] * 0)
    ring.AddPoint(geotf[0] + geotf[1] * xsize, geotf[3] + geotf[5] * ysize)
    ring.AddPoint(geotf[0] + geotf[1] * 0, geotf[3] + geotf[5] * ysize)
    ring.AddPoint(geotf[0] + geotf[1] * 0, geotf[3] + geotf[5] * 0)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly.AssignSpatialReference(srs)
    poly.TransformTo(layer.GetSpatialRef())
    env = poly.GetEnvelope()
    layer.SetSpatialFilterRect(env[0], env[2], env[1], env[3])
    # Rasterize
    driver = gdal.GetDriverByName('MEM')
    dataset = driver.Create('', xsize, ysize, 1, gdal.GDT_Byte)
    dataset.SetProjection(srs.ExportToWkt())
    dataset.SetGeoTransform(geotf)
    gdal.RasterizeLayer(dataset, [1], layer, burn_values=[1],
                        options=['ALL_TOUCHED=TRUE'])
    landmask = dataset.GetRasterBand(1).ReadAsArray().astype('bool')
    return landmask


def make_datagroup(stitching_groups):
    """
    """
    # Get SAFE name for each granule
    safe_names = []
    for granule_paths in stitching_groups.values():
        safe_paths = [os.path.dirname(os.path.dirname(path)) for path in granule_paths]
        safe_names.extend([os.path.basename(path) for path in safe_paths])
    # If multiple SAFE names (single tile products), we mix first and last one in lexical order.
    safe_names = list(set(safe_names))
    if len(safe_names) > 1:
        safe_names.sort()
        safe_split1 = safe_names[0].split('_')
        safe_split2 = safe_names[-1].split('_')
        safe_split = []
        for f1, f2 in zip(safe_split1, safe_split2):
            if f1 == f2:
                safe_split.append(f1)
            else:
                safe_split.append('{}-{}'.format(f1, f2))
        safe_name = '_'.join(safe_split)
    else:
        safe_name = safe_names[0]
    datagroup = os.path.splitext(safe_name)[0]
    return datagroup


def set_contrast(vmin, vmax, groups, overview_index=2, landmaskpath=None,
                 slope_threshold=-40., debug_fig_dir=None, atmos_correction=0,
                 atmos_lut_path=None):
    """
    """
    # Get bands for each projection
    bandnames = ['B04', 'B03', 'B02']
    bands = {}
    for proj, urls in groups.iteritems():
        # Mapper
        print '    {} : open mapper with overview {}'.format(proj, overview_index)
        t0 = datetime.utcnow()
        mapper = SAFEMSIL1CStitchedFile(urls, native_resolution='10m',
                                        overview_index=overview_index, tight=True,
                                        atmos_lut_path=atmos_lut_path)
        mapper.open()
        print '        {}'.format(datetime.utcnow() - t0)
        _bands = {}
        spacing = [abs(mapper.read_global_attribute('ydim')),
                   abs(mapper.read_global_attribute('xdim'))]
        shape = (mapper.get_dimsize('y'), mapper.get_dimsize('x'))
        _bands['qvalue'] = mapper.read_global_attribute('quantification_value')
        # Masks
        print '    {} : get landmask'.format(proj)
        t0 = datetime.utcnow()
        if landmaskpath is not None and os.path.exists(landmaskpath):
            landmask = compute_landmask(mapper, landmaskpath, downsampling=1)
            dil_len = 500.
            dil_shp = [2 * np.round(dil_len / sp).astype('int') + 1 for sp in spacing]
            dil_kern = np.ones(dil_shp, dtype='bool')
            landmask = binary_dilation(landmask, structure=dil_kern)
            _bands['landmask'] = landmask
        else:
            print 'WARNING : no land mask path was given.'
            _bands['landmask'] = np.zeros(shape, dtype='bool')
        print '        {}'.format(datetime.utcnow() - t0)
        print '    {} : get cloudmask'.format(proj)
        t0 = datetime.utcnow()
        if 'cloud_mask' in mapper.get_fieldnames():
            cloudmask = mapper.read_values('cloud_mask')
            dil_len = 500.
            dil_shp = [2 * np.round(dil_len / sp).astype('int') + 1 for sp in spacing]
            dil_kern = np.ones(dil_shp, dtype='bool')
            cloudmask = binary_dilation(cloudmask, structure=dil_kern)
            _bands['cloudmask'] = cloudmask
        else:
            print 'WARNING : no cloud mask in product.'
            _bands['cloudmask'] = np.zeros(shape, dtype='bool')
        print '        {}'.format(datetime.utcnow() - t0)
        # Bands
        for bandname in bandnames:
            fieldname = '{}_digital_number'.format(bandname)
            print '    {} : read {}'.format(proj, fieldname)
            t0 = datetime.utcnow()
            band_dn = mapper.read_values(fieldname)
            if isinstance(band_dn, np.ma.MaskedArray):
                _bands[bandname] = band_dn
            else:
                _bands[bandname] = np.ma.MaskedArray(band_dn)
            print '        {}'.format(datetime.utcnow() - t0)
            if atmos_correction != 0:
                fieldname = '{}_atmospheric_digital_number'.format(bandname)
                print '    {} : read {}'.format(proj, fieldname)
                t0 = datetime.utcnow()
                band_atm_dn = mapper.read_values(fieldname)
                if isinstance(band_atm_dn, np.ma.MaskedArray):
                    _bands['atm_{}'.format(bandname)] = band_atm_dn
                else:
                    _bands['atm_{}'.format(bandname)] = np.ma.MaskedArray(band_atm_dn)
                print '        {}'.format(datetime.utcnow() - t0)
        bands[proj] = _bands
        mapper.close()

    # Get valid values for each band
    print '    Concatenate valid values'
    t0 = datetime.utcnow()
    values = {}
    for bandname in bandnames:
        _data = [d[bandname].data for d in bands.values()]
        _mask = [d[bandname].mask | d['landmask'] | d['cloudmask'] for d in bands.values()]
        values[bandname] = np.concatenate([_d[~_m] for _d, _m in zip(_data, _mask)])
    if atmos_correction != 0:
        for bandname in bandnames:
            _bandname = 'atm_{}'.format(bandname)
            _data = [d[_bandname].data for d in bands.values()]
            _mask = [d[_bandname].mask | d['landmask'] | d['cloudmask'] for d in bands.values()]
            values[_bandname] = np.concatenate([_d[~_m] for _d, _m in zip(_data, _mask)])
    print '        {}'.format(datetime.utcnow() - t0)
    qvalue = [b['qvalue'] for b in bands.values()]
    if min(qvalue) != max(qvalue):
        raise Exception('Groups do not share the same quantification value ?')
    qvalue = qvalue[0]

    # Set vmin
    print '    Set vmin'
    t0 = datetime.utcnow()
    if atmos_correction == 0:
        for iband, bandname in enumerate(bandnames):
            if vmin[iband] is not None:
                continue
            vmin[iband] = np.percentile(values[bandname], 0.5) / qvalue
    elif atmos_correction == 1:
        for iband, bandname in enumerate(bandnames):
            if vmin[iband] is not None:
                continue
            vmin[iband] = values['atm_{}'.format(bandname)].mean() / qvalue
    elif atmos_correction == 2:
        _vmin, _vmin_atm = [], []
        for iband, bandname in enumerate(bandnames):
            _vmin.append(np.percentile(values[bandname], 0.5) / qvalue)
            _vmin_atm.append(values['atm_{}'.format(bandname)].mean() / qvalue)
        shifts = [_v - _v_atm for _v, _v_atm in zip(_vmin, _vmin_atm)]
        shift = min(shifts)
        for iband, bandname in enumerate(bandnames):
            if vmin[iband] is not None:
                continue
            vmin[iband] = _vmin_atm[iband] + shift
    print '        {}'.format(datetime.utcnow() - t0)

    # Set vmax
    print '    Set vmax'
    t0 = datetime.utcnow()
    if atmos_correction == 0:
        for iband, bandname in enumerate(bandnames):
            if vmax[iband] is not None:
                continue
            if iband == 0: # red
                bins = np.linspace(0, .5 * qvalue, num=51)
                hist, edges = np.histogram(values[bandname], bins=bins)
                edges = edges / qvalue
                dbin = edges[1] - edges[0]
                hist = hist.astype('float') / hist.sum() / dbin
                slope = (hist[1:] - hist[:-1]) / dbin
                slope_ctrs = edges[1:-1]
                ind = np.where(slope <= slope_threshold)[0].max()
                dy = slope[ind + 1] - slope[ind]
                dx = slope_ctrs[ind + 1] - slope_ctrs[ind]
                vmax[iband] = slope_ctrs[ind] + (slope_threshold - slope[ind]) / dy * dx
            elif iband == 1 and atmos_correction == 0: # green
                vmax[iband] = vmin[iband] + (vmax[0] - vmin[0]) * 0.887
            elif iband == 2 and atmos_correction == 0: # blue
                vmax[iband] = vmin[iband] + (vmax[0] - vmin[0]) * 0.907
    elif atmos_correction > 0:
        _vmax = []
        for iband, bandname in enumerate(bandnames):
            bins = np.linspace(0, .5 * qvalue, num=51)
            hist, edges = np.histogram(values[bandname], bins=bins)
            edges = edges / qvalue
            dbin = edges[1] - edges[0]
            hist = hist.astype('float') / hist.sum() / dbin
            slope = (hist[1:] - hist[:-1]) / dbin
            slope_ctrs = edges[1:-1]
            ind = np.where(slope <= slope_threshold)[0].max()
            dy = slope[ind + 1] - slope[ind]
            dx = slope_ctrs[ind + 1] - slope_ctrs[ind]
            _vmax.append(slope_ctrs[ind] + (slope_threshold - slope[ind]) / dy * dx)
        vrans = [vma - vmi for vma, vmi in zip(_vmax, vmin)]
        vran = max(vrans)
        for iband, bandname in enumerate(bandnames):
            if vmax[iband] is not None:
                continue
            vmax[iband] = vmin[iband] + vran
    print '        {}'.format(datetime.utcnow() - t0)

    # Make figures for debug
    if debug_fig_dir is not None and os.path.exists(debug_fig_dir):
        print '    Make debug figures'
        t0 = datetime.utcnow()
        import matplotlib.pyplot as plt
        datagroup = make_datagroup(groups)
        # Show bands and masks
        for proj, _bands in bands.iteritems():
            fig = plt.figure(figsize=(12, 9))
            for iband, bandname in enumerate(bandnames):
                plt.subplot(2, 2, iband + 1)
                plt.imshow(_bands[bandname], vmin=vmin[iband] * qvalue,
                           vmax=vmax[iband] * qvalue)
                plt.colorbar(shrink=0.5)
                plt.title(bandname)
                plt.xticks([])
                plt.yticks([])
            plt.subplot(2, 2, 4)
            mask = _bands['landmask'] | _bands['cloudmask']
            plt.imshow(mask, vmin=0, vmax=1)
            plt.colorbar(shrink=0.5)
            plt.title('land/cloud mask')
            plt.xticks([])
            plt.yticks([])
            plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975,
                                wspace=0.05, hspace=0.05)
            figpath = os.path.join(debug_fig_dir, '{}-{}.png'.format(datagroup, proj[5:]))
            plt.savefig(figpath)
            plt.close(fig)
        if atmos_correction != 0:
            for proj, _bands in bands.iteritems():
                fig = plt.figure(figsize=(12, 9))
                for iband, bandname in enumerate(bandnames):
                    plt.subplot(2, 2, iband + 1)
                    _bandname = 'atm_{}'.format(bandname)
                    plt.imshow(_bands[_bandname])
                    plt.colorbar(shrink=0.5)
                    plt.title(_bandname)
                    plt.xticks([])
                    plt.yticks([])
                plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975,
                                    wspace=0.05, hspace=0.05)
                figpath = os.path.join(debug_fig_dir, '{}-{}-atm.png'.format(datagroup, proj[5:]))
                plt.savefig(figpath)
                plt.close(fig)
        # Show histograms
        colors = ['red', 'green', 'blue']
        fig = plt.figure(figsize=(12, 9))
        plt.subplot(2, 1, 1)
        for iband, bandname in enumerate(bandnames):
            bins = np.linspace(0, .5 * qvalue, num=501)
            hist, edges = np.histogram(values[bandname], bins=bins)
            edges = edges / qvalue
            dbin = edges[1] - edges[0]
            hist = hist.astype('float') / hist.sum() / dbin
            hist_ctrs = edges[:-1] + dbin / 2.
            plt.plot(hist_ctrs, hist, color=colors[iband])
        ylim = plt.ylim()
        vmintxt, vmaxtxt, vrantxt = 'vmin=', 'vmax=', 'vmax-vmin='
        for iband, bandname in enumerate(bandnames):
            plt.plot([vmin[iband]] * 2, ylim, '--', color=colors[iband])
            plt.plot([vmax[iband]] * 2, ylim, '--', color=colors[iband])
            vmintxt += '{:.3f} '.format(vmin[iband])
            vmaxtxt += '{:.3f} '.format(vmax[iband])
            vrantxt += '{:.3f} '.format(vmax[iband] - vmin[iband])
        plt.xlabel('TOA reflectance')
        plt.ylabel('Histogram')
        plt.suptitle(' / '.join([vmintxt, vmaxtxt, vrantxt]))
        plt.subplot(2, 1, 2)
        for iband, bandname in enumerate(bandnames):
            bins = np.linspace(0, .5 * qvalue, num=51)
            hist, edges = np.histogram(values[bandname], bins=bins)
            edges = edges / qvalue
            dbin = edges[1] - edges[0]
            hist = hist.astype('float') / hist.sum() / dbin
            slope = (hist[1:] - hist[:-1]) / dbin
            slope_ctrs = edges[1:-1]
            plt.plot(slope_ctrs, slope, '+-', color=colors[iband])
        #plt.ylim(-50, 50)
        plt.ylim(-200, 10)
        ylim = plt.ylim()
        for iband, bandname in enumerate(bandnames):
            #plt.plot([vmin[iband]] * 2, ylim, '--', color=colors[iband])
            plt.plot([vmax[iband]] * 2, ylim, '--', color=colors[iband])
        xlim = plt.xlim()
        plt.plot(xlim, [slope_threshold] * 2, '-r')
        plt.xlabel('TOA reflectance')
        plt.ylabel('Histogram slope')
        figpath = os.path.join(debug_fig_dir, '{}-{}.png'.format(datagroup, 'histo'))
        plt.savefig(figpath)
        plt.close(fig)
        print '        {}'.format(datetime.utcnow() - t0)

    return vmin, vmax


def sentinel2_rgb(infile, outdir,
                  # For output resolution
                  overview_index=None, downsampling=2,
                  # For manual contrast
                  vmin=[None, None, None], vmax=[None, None, None],
                  # For auto contrast
                  contrast_overview_index=2, landmaskpath=None,
                  slope_threshold=-40., debug_fig_dir=None,
                  atmos_correction=0, atmos_lut_path=None,
                  vmax_factor=None,
                  # For output type
                  write_netcdf=False):
    """
    """
    # Identify stitching groups
    print 'Identify stitching group(s)'
    groups = safemsil1c_stitching_groups(infile)
    projs = groups.keys()
    for proj, urls in groups.iteritems():
        print '    {} : {} granule(s)'.format(proj, len(urls))
    datagroup = make_datagroup(groups)

    # Set contrast
    if None in vmin or None in vmax:
        print 'Set contrast'
        vmin, vmax = set_contrast(list(vmin), list(vmax), groups,
                                  overview_index=contrast_overview_index,
                                  landmaskpath=landmaskpath,
                                  slope_threshold=slope_threshold,
                                  debug_fig_dir=debug_fig_dir,
                                  atmos_correction=atmos_correction,
                                  atmos_lut_path=atmos_lut_path)
    print 'vmin = {:0.4f} / {:0.4f} / {:0.4f}'.format(*vmin)
    print 'vmax = {:0.4f} / {:0.4f} / {:0.4f}'.format(*vmax)
    if vmax_factor is not None:
        print 'Apply vmax_factor={}'.format(vmax_factor)
        _vmax = []
        for vmi, vma in zip(vmin, vmax):
            _vmax.append(vmi + vmax_factor * (vma - vmi))
        vmax = _vmax
        print 'new vmax = {:0.4f} / {:0.4f} / {:0.4f}'.format(*vmax)

    # Build geotiff or netcdf
    print 'Build geotiff or netcdf'
    bandnames = ['B04', 'B03', 'B02']
    for proj in projs:
        # Open stitched mapper
        print '    {} : open mapper with overview {}'.format(proj, overview_index)
        t0 = datetime.utcnow()
        mapper = SAFEMSIL1CStitchedFile(groups[proj],
                                        native_resolution='10m',
                                        overview_index=overview_index,
                                        tight=True)
        mapper.open()
        print '        {}'.format(datetime.utcnow() - t0)

        # Construct bands
        bands = []
        qvalue = mapper.read_global_attribute('quantification_value')
        for iband, bandname in enumerate(bandnames):
            fieldname = '{}_digital_number'.format(bandname)
            print '    {} : read {}'.format(proj, fieldname)
            t0 = datetime.utcnow()
            band = mapper.read_values(fieldname)
            print '        {}'.format(datetime.utcnow() - t0)
            if downsampling != 1:
                print '    {} : downsample by {}'.format(proj, downsampling)
                t0 = datetime.utcnow()
                shp = list(band.shape)
                shp[0] -= np.mod(shp[0], downsampling)
                shp[1] -= np.mod(shp[1], downsampling)
                sli = [slice(0, shp[0]), slice(0, shp[1])]
                rshp = (shp[0] / downsampling, downsampling,
                        shp[1] / downsampling, downsampling)
                if not np.ma.is_masked(band):
                    mask = np.ma.nomask
                else:
                    mask = band[sli].mask.reshape(rshp).\
                           sum(axis=3, dtype='uint8').\
                           sum(axis=1, dtype='uint8') > 0
                band = np.ma.MaskedArray(band[sli].data.reshape(rshp).\
                                         mean(axis=3, dtype='uint16').\
                                         mean(axis=1, dtype='uint16'),
                                         mask=mask)
                del mask
                print '        {}'.format(datetime.utcnow() - t0)
            print '    {} : bytescale in [{}, {}]'.format(proj, vmin[iband], vmax[iband])
            t0 = datetime.utcnow()
            vmin_dn = np.round(vmin[iband] * qvalue)
            vmax_dn = np.round(vmax[iband] * qvalue)
            byte = bytescale(band.data, cmin=vmin_dn, cmax=vmax_dn,
                             low=0, high=254)
            if band.mask is not np.ma.nomask:
                byte[band.mask] = 255
            del band
            scale = (vmax[iband] - vmin[iband]) / 254.
            offset = vmin[iband]
            description = '{} TOA reflectance'.format(bandname)
            bands.append({'array': byte,
                          'scale': scale,
                          'offset': offset,
                          'description': description,
                          'unittype': '',
                          'nodatavalue': 255,
                          'parameter_range': [vmin[iband], vmax[iband]]})
            print '        {}'.format(datetime.utcnow() - t0)

        # Make sure nodata are at the same locations in all bands
        mask = np.any([band['array'] == 255 for band in bands], axis=0)
        for band in bands:
            band['array'][mask] = 255

        # Construct metadata and geolocation
        print '    {} : construct metadata and geolocation'.format(proj)
        t0 = datetime.utcnow()
        cs_code = mapper.read_global_attribute('horizontal_cs_code')
        epsg_num = cs_code.lower().lstrip('epsg:')
        dataname = '{}-{}'.format(datagroup, epsg_num)
        start_time = mapper.get_start_time()
        end_time = mapper.get_end_time()
        (dtime, time_range) = stfmt.format_time_and_range(start_time, end_time,
                                                          units='ms')
        sensor_pass = mapper.read_global_attribute('sensing_orbit_direction').lower()
        metadata = {}
        metadata['product_name'] = 'Sentinel-2_RGB'
        metadata['name'] = dataname
        metadata['datetime'] = dtime
        metadata['time_range'] = time_range
        metadata['source_URI'] = infile
        metadata['source_provider'] = 'ESA'
        metadata['processing_center'] = 'OceanDataLab'
        metadata['conversion_software'] = 'Syntool'
        metadata['conversion_version'] = '0.0.0'
        metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
        metadata['parameter'] = ['B04 TOA reflectance',
                                 'B03 TOA reflectance',
                                 'B02 TOA reflectance']
        metadata['type'] = 'remote sensing'
        metadata['sensor_type'] = 'multi-spectral imager'
        metadata['sensor_name'] = 'MSI'
        metadata['sensor_platform'] = 'Sentinel-2'
        metadata['sensor_pass'] = sensor_pass
        metadata['datagroup'] = datagroup
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg_num))
        ulx = mapper.read_global_attribute('ulx')
        dx = mapper.read_global_attribute('xdim') * downsampling
        uly = mapper.read_global_attribute('uly')
        dy = mapper.read_global_attribute('ydim') * downsampling
        geolocation = {}
        geolocation['projection'] = srs.ExportToWkt()
        geolocation['geotransform'] = [ulx, dx, 0, uly, 0, dy]
        print '        {}'.format(datetime.utcnow() - t0)

        # Write geotiff or netcdf
        mapper.close()
        if write_netcdf == False:
            print '    {} : write geotiff'.format(proj)
            t0 = datetime.utcnow()
            tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
            stfmt.write_geotiff(tifffile, metadata, geolocation, bands)
            print '        {}'.format(datetime.utcnow() - t0)
        elif write_netcdf == True:
            print '    {} : write geotiff'.format(proj)
            t0 = datetime.utcnow()
            ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
            bands[0]['name'] = 'B04_TOA_reflectance'
            bands[1]['name'] = 'B03_TOA_reflectance'
            bands[2]['name'] = 'B02_TOA_reflectance'
            resolution = min([abs(dy), abs(dx)])
            metadata['spatial_resolution'] = resolution
            dgcps_meter = 25000.
            dgcps = (np.round(dgcps_meter / resolution).astype('int'), ) * 2
            stfmt.write_netcdf(ncfile, metadata, geolocation, bands, 'grid_proj',
                               dgcps=dgcps)
            print '        {}'.format(datetime.utcnow() - t0)
