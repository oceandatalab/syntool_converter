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
from syntool_converter.utils.sst import denoise_sst
import os
import pyproj
import re
import sys

# def resampling_test():
#     """ """
#     import gdal
#     import matplotlib.pyplot as plt
#     mrtswathsstfname = '/mnt/data/modis/resampling_tests/A2011338122500.L2_LAC_SST_sst.tif'
#     dset = gdal.Open(mrtswathsstfname)
#     sst = dset.GetRasterBand(1).ReadAsArray()
#     gdalproj = dset.GetProjection()
#     gdalgeot = dset.GetGeoTransform()
#     dset = None
#     mask = (sst == 0) | (sst == -32767)
#     sst = sst * 0.0049999999 + 273.5
#     #mrtswathcloudmaskfname = '/mnt/data/modis/resampling_tests/MYD35_L2.A2011338.1225.006.2014086060029_Cloud_Mask_b0.tif'
#     mrtswathcloudmaskfname = '/mnt/data/modis/resampling_tests/MYD35_L2.A2011338.1225.005.2011340001234_Cloud_Mask_b0.tif'
#     dset = gdal.Open(mrtswathcloudmaskfname)
#     cloudmask = dset.GetRasterBand(1).ReadAsArray()
#     dset = None
#     valid = np.where(mask == False)
#     cloudy = (np.bitwise_and(cloudmask[valid], 2) == 0) & \
#              (np.bitwise_and(cloudmask[valid], 4) == 0)
#     land = np.bitwise_and(cloudmask[valid], 128) == 128
#     mask[valid] = cloudy | land
#     #import pdb ; pdb.set_trace()
#     metadata = {}
#     geolocation = {}
#     geolocation['projection'] = gdalproj
#     geolocation['geotransform'] = gdalgeot
#     band = []
#     vmin = 283.511330872
#     vmax = 298.79493927
#     indndv = np.where(mask == True)
#     offset, scale = vmin, (vmax-vmin)/254.
#     np.clip(sst, vmin, vmax, out=sst)
#     array = np.round((sst - offset) / scale).astype('uint8')
#     array[indndv] = 255
#     colortable = stfmt.format_colortable('matplotlib_jet',
#                                          vmax=vmax, vmax_pal=vmax,
#                                          vmin=vmin, vmin_pal=vmin)
#     band.append({'array':array, 'scale':scale, 'offset':offset,
#                  'description':'sea surface temperature', 'unittype':'K',
#                  'nodatavalue':255, 'parameter_range':[vmin, vmax],
#                  'colortable':colortable})
#     tifffile = '/mnt/data/modis/resampling_tests/sst_mrtswath.tiff'
#     stfmt.write_geotiff(tifffile, metadata, geolocation, band)


# swath size : ~2030km x ~2357.05km
# 250km ngcps=(9, 10)
# 200km ngcps=(11, 13)
# 150km ngcps=(15, 17)
# 100km ngcps=(21, 25)
# 075km ngcps=(28, 32)
# 050km ngcps=(42, 48)
# 025km ngcps=(82, 95)
# modis_sst(infileid='A2011338122500', outdir='/mnt/data/syntool_inputs/synergy/',
#           download_dir='/mnt/data/modis')

def modis_sst(infileid, outdir, download_dir='/tmp',
              vmin=None, vmax=None, contrast='relative',
              ngcps=(21, 25), resample_radius=5000., resample_sigma=2500.,
              denoise_kernel='boxcar', denoise_width=20, open_iterations=1,
              nprocs=1, pngkml=False, write_netcdf=False, file_range=None):
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

    # modissstfname = '/mnt/data/sst/modis/MYD021KM.A2011338.1225/A2011338122500.L2_LAC_SST'
    # modis02fname = '/mnt/data/sst/modis/MYD021KM.A2011338.1225/MYD021KM.A2011338.1225.005.2011339235825.hdf'
    # modis03fname = '/mnt/data/sst/modis/MYD021KM.A2011338.1225/MYD03.A2011338.1225.005.2011339233301.hdf'
    # modis35l2fname = '/mnt/data/sst/modis/MYD021KM.A2011338.1225/MYD35_L2.A2011338.1225.005.2011340001234.hdf'
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
        raise Exception('Input for modis_sst is an ID '
                        '(e.g. A2011338122500 or T2014143234500)')
    platform = infileid[0]
    date = datetime.strptime(infileid[1:], '%Y%j%H%M%S')
    modissstid = {'A':'MODISAL2SST', 'T':'MODISTL2SST'}[platform]
    modissstfname = modis.search_and_download(modissstid, date, download_dir)
    modis02id = {'A':'MYD021KM', 'T':'MOD021KM'}[platform]
    modis02fname = modis.search_and_download(modis02id, date, download_dir)
    modis03id = {'A':'MYD03', 'T':'MOD03'}[platform]
    modis03fname = modis.search_and_download(modis03id, date, download_dir)
    modis35l2id = {'A':'MYD35_L2', 'T':'MOD35_L2'}[platform]
    modis35l2fname = modis.search_and_download(modis35l2id, date, download_dir)

    # Read/Process data
    print 'Read/Process data'
    # Read from SST file
    modissstfile = modis.MODISL2File(modissstfname)
    # lon = modissstfile.read_lon()
    # lat = modissstfile.read_lat()
    sst = modissstfile.read_sst() + 273.15
    attrs = modissstfile.read_attributes()
    modissstfile.close()
    # Read from radiances file
    modis02file = modis.MODIS02File(modis02fname)
    rad11 = modis02file.read_radiance(31)
    modis02file.close()
    bt11 = modis.modis_bright(rad11, 31, 1)
    # Read from geolocation file
    modis03file = modis.MODIS03File(modis03fname)
    lon = modis03file.read_lon()
    lat = modis03file.read_lat()
    modis03file.close()
    # Read from cloud mask file
    modis35l2file = modis.MODIS35L2File(modis35l2fname)
    cloudmask = modis35l2file.read_cloudmask(byte=0)
    modis35l2file.close()
    cloudy = (np.bitwise_and(cloudmask, 2) == 0) & \
             (np.bitwise_and(cloudmask, 4) == 0)
    land = np.bitwise_and(cloudmask, 128) == 128 # Desert or Land
    # land = (np.bitwise_and(cloudmask, 128) == 128) | \
    #        (np.bitwise_and(cloudmask, 64) == 64) # Desert or Land or Coastal
    if listbox is not None:
        mask_box = np.zeros(np.shape(sst))
        for i in range(np.shape(listbox)[0]):
            index_in = np.where((lon >= listbox[i][0]) & (lat >= listbox[i][1])
                             & (lon <= listbox[i][2]) & (lat <= listbox[i][3]))
            mask_box[index_in] = 1
        mask = cloudy | land | ma.getmaskarray(sst) | ma.getmaskarray(bt11)| (mask_box == 0)
    else:
        mask = cloudy | land | ma.getmaskarray(sst) | ma.getmaskarray(bt11)
    if mask.all():
        print 'No data'
        sys.exit(0)
    # GCPs for resampling and geotiff georeference
    scansize = 10
    dtime0 = datetime.utcnow()
    gcps = resample.get_gcps_from_bowtie(lon, lat, scansize, ngcps=ngcps)
    #gcps = resample.get_gcps_from_bowtie_old(lon, lat, scansize, ngcps=ngcps)
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

    # Resample with pyresample in lon/lat space
    # rsplin, rsppix = np.mgrid[0:rspysize, 0:rspxsize] + 0.5
    # rsplon, rsplat = resample.get_points_from_gcps(gcplon, gcplat, gcppixel,
    #                                                gcpline, rspxsize, rspysize,
    #                                                0, rsppix, rsplin, nprocs=nprocs)
    # # Test resample grid
    # import matplotlib.pyplot as plt
    # plt.plot(lon.flatten(), lat.flatten(), '+b')
    # plt.plot(rsplon.flatten(), rsplat.flatten(), '+g')
    # plt.plot(gcplon.flatten(), gcplat.flatten(), 'xr')
    # plt.show()
    # import pdb ; pdb.set_trace()
    # # \Test resample grid
    # # Test radius / sigma
    # resample_radius = 5000.
    # resample_sigma = 2500.
    # sst.mask = False
    # #sst.mask = sst.mask | (sst.data < 273.15+5) | (sst.data > 273.15+30)
    # rspsst = resample.resample_gauss(lon, lat, sst, rsplon, rsplat,
    #                                  resample_radius, resample_sigma,
    #                                  nprocs=nprocs, show=True)
    # import pdb ; pdb.set_trace()
    # # \Test radius / sigma
    # valid = np.where(mask == False)
    # rspsst = resample.resample_gauss(lon[valid], lat[valid], sst[valid],
    #                                  rsplon, rsplat,
    #                                  resample_radius, resample_sigma,
    #                                  fill_value=None, nprocs=nprocs,
    #                                  show=False)
    # rspbt11 = resample.resample_gauss(lon[valid], lat[valid], bt11[valid],
    #                                   rsplon, rsplat,
    #                                   resample_radius, resample_sigma,
    #                                   fill_value=None, nprocs=nprocs,
    #                                   show=False)
    # rspmask = resample.resample_nearest(lon, lat, mask,
    #                                     rsplon, rsplat,
    #                                     resample_radius,
    #                                     fill_value=True, nprocs=nprocs,
    #                                     show=False)
    # rspmask = rspmask | ma.getmaskarray(rspsst) | ma.getmaskarray(rspbt11)

    # Denoise sst and open mask
    rspsst.mask = rspmask
    rspbt11.mask = rspmask
    finalsst = denoise_sst(rspsst, rspbt11, kernel=denoise_kernel,
                           width=denoise_width, show=False)
    #finalsst = rspsst
    finalmask = ~binary_opening(~rspmask, structure=np.ones((3, 3)),
                                iterations=open_iterations)
    #finalmask = rspmask
    finalsst.mask = finalmask

    # Contrast
    if vmin == None:
        if contrast == 'relative':
            vmin = np.percentile(finalsst.compressed(), 0.5)
        #elif contrast == 'agulhas':
        #    dayofyear = float(attrs['start_time'].timetuple().tm_yday)
        #    vmin = 273.15 + 2. * np.cos((dayofyear - 45.) * 2. * np.pi / 365.) + 20. - 9.
        #    #par = [277.94999694824219, 42, 2.5500030517578125, -219]
        #    par = [278.09999084472656, 0.62831853071795862,
        #           2.4000091552734375, 0.1570796326794896]
        #    vmin = par[0] + par[2] * np.cos(par[3] * dayofyear - par[1])
        #if a specific txt file is provided for the range
        elif dic is not None:
            dayofyear = float(attrs['start_time'].timetuple().tm_yday)
            # Read a txt file which contains three columns: yearday,vmin,vmax
            extrema = dic.get(dayofyear, dic[min(dic.keys(),
                           key=lambda k:abs(k - dayofyear))])
            vmin = extrema[0]
        else:
            raise Exception('Unknown contrast : {}'.format(contrast))
    if vmax == None:
        if contrast == 'relative':
            vmax = np.percentile(finalsst.compressed(), 99.5)
        #elif contrast == 'agulhas':
        #    dayofyear = float(attrs['start_time'].timetuple().tm_yday)
        #    vmax = 273.15 + 2. * np.cos((dayofyear - 45.) * 2. * np.pi / 365.) + 20. + 4.
        #    #par = [300.59999084472656, 21, 2.8499908447265625, -191]
        #    par = [300.59999084472656, 0.29919930034188508,
        #           2.8499908447265625, 0.14959965017094254]
        #    vmax = par[0] + par[2] * np.cos(par[3] * dayofyear - par[1])
        #if a specific text file is provided for the range
        elif dic is not None:
            dayofyear = float(attrs['start_time'].timetuple().tm_yday)
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
    (dtime, time_range) = stfmt.format_time_and_range(attrs['start_time'],
                                                      attrs['stop_time'],
                                                      units='ms')
    metadata['product_name'] = 'SST_MODIS_denoised'
    if contrast == 'relative':
        metadata['name'] = os.path.splitext(os.path.basename(modissstfname))[0]
    else:
        metadata['name'] = '{}_{}'.format(os.path.splitext(os.path.basename(modissstfname))[0], contrast)
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = [modissstfname, modis02fname,
                              modis03fname, modis35l2fname]
    metadata['source_provider'] = 'NASA'
    metadata['processing_center'] = 'OceanDataLab'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = 'sea surface temperature'
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
        # # e.g. 999.763079208 999.763084628 999.763082543
        # print ydists.min(), ydists.max(), ydists.mean()
        # # e.g. 1006.4149472 1008.60679776 1007.5888004
        metadata['spatial_resolution'] = 1000.
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'swath',
                           ngcps=gcplon.shape)

