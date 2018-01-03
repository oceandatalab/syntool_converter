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

from cerbere.mapper.ecmwf0125ncfile import ECMWF0125NCFile
import numpy as np
from netCDF4 import num2date, Dataset
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import os


# def make_land_mask():
#     """
#     """
#     from cerbere.mapper.ncfile import NCFile
#     lmfile = NCFile('/local/home/data/ancillary/land_mask/landmask_025.nc')
#     mask = lmfile.read_values('mask')[0, :, :]
#     # resample #1
#     # from pyresample import geometry, kd_tree
#     # lon = lmfile.read_values('lon')
#     # lat = lmfile.read_values('lat')
#     # newlon = np.arange(2880)*0.125 - 180.
#     # newlat = np.arange(1409)*0.125 - 88.
#     # lonG = np.tile(lon.reshape((1, -1)), (lat.size, 1))
#     # latG = np.tile(lat.reshape((-1, 1)), (1, lon.size))
#     # grid = geometry.GridDefinition(lonG, latG)
#     # newlonG = np.tile(newlon.reshape((1, -1)), (newlat.size, 1))
#     # newlatG = np.tile(newlat.reshape((-1, 1)), (1, newlon.size))
#     # newgrid = geometry.GridDefinition(newlonG, newlatG)
#     # newmask = kd_tree.resample_nearest(grid, mask, newgrid, 100000.)
#     # resample #2
#     mask = mask[8:713, :]
#     newmask = np.repeat(mask, 2, axis=0)
#     newmask = np.repeat(newmask, 2, axis=1)
#     newmask = newmask[0:-1, :]
#     # save
#     fname = os.path.join(os.path.dirname(__file__), 'ecmwf_0125_land_mask.nc')
#     rootgrp = Dataset(fname, 'w', format='NETCDF4')
#     lat = rootgrp.createDimension('lat', 1409)
#     lon = rootgrp.createDimension('lon', 2880)
#     mask = rootgrp.createVariable('mask', 'u1', ('lat', 'lon'), zlib=True)
#     mask[:, :] = np.array(newmask, dtype='uint8')
#     rootgrp.close()


def get_land_mask():
    """
    """
    fname = os.path.join(os.path.dirname(__file__), 'ecmwf_0125_land_mask.nc')
    rootgrp = Dataset(fname, 'r', format='NETCDF4')
    land_mask = rootgrp.variables['mask'][:, :]
    rootgrp.close()
    return land_mask


def ecmwf_model_wind(infile, outdir, max_forecast_hours=None,
                     vmin=0., vmax=25.4, vmin_pal=0., vmax_pal=50*0.514,
                     write_netcdf=False):
    """
    """
    # Read/Process data
    windfield = ECMWF0125NCFile(infile)
    u10 = windfield.read_values('u10m')[0, ::-1, :]
    v10 = windfield.read_values('v10m')[0, ::-1, :]
    lon = windfield.read_values('lon')
    dlon = lon[1]-lon[0]
    lat = windfield.read_values('lat')[::-1]
    dlat = lat[1]-lat[0]
    land_mask = get_land_mask()[::-1, :]
    # Replicate -180 deg at 180 deg for gdal_warp
    # dim = u10.shape
    # u10 = np.hstack((u10, u10[:, 0].reshape((dim[0], 1))))
    # v10 = np.hstack((v10, v10[:, 0].reshape((dim[0], 1))))
    # lon = np.hstack((lon, lon[0]+360.))
    # land_mask = np.hstack((land_mask, land_mask[:, 0].reshape((dim[0], 1))))
    # /Replicate -180 deg at 180 deg for gdal_warp
    dtime = windfield.read_values('time')[0]
    dtime_units = windfield.read_field('time').units
    dtime = num2date(dtime, dtime_units)
    rundtime = windfield.read_global_attribute('run_time')
    rundtime = datetime.strptime(rundtime, '%Y-%m-%dT%H:%M:%SZ')
    if max_forecast_hours is not None:
        forecast_hours = (dtime - rundtime).total_seconds() / 3600.
        if forecast_hours > max_forecast_hours:
            raise Exception('Exceeds max_forecast_hours.')
    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    metadata['product_name'] = 'ECMWF_model_wind'
    #metadata['name'] = os.path.splitext(os.path.basename(infile))[0]
    metadata['name'] = 'ECMWF_'+dtime.strftime('%Y%m%dT%HZ')
    metadata['datetime'] = stfmt.format_time(dtime)
    metadata['time_range'] = ['-90m', '+90m']
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'ECMWF'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    #metadata['parameter'] = ['zonal wind speed', 'meridional wind speed']
    metadata['parameter'] = ['wind speed', 'wind direction']
    metadata['type'] = 'model'
    metadata['model_longitude_resolution'] = 0.125
    metadata['model_latitude_resolution'] = 0.125
    metadata['model_analysis_datetime'] = stfmt.format_time(rundtime)
    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['geotransform'] = [lon[0]-dlon/2., dlon, 0,
                                   lat[0]-dlat/2., 0, dlat]
    # band = []
    # scale = 0.2
    # offset = -25.4
    # windspeed = np.sqrt(u10**2 + v10**2)
    # winddirection = np.arctan2(v10, u10)
    # np.clip(windspeed, 0, abs(offset), out=windspeed)
    # u10 = np.cos(winddirection)*windspeed
    # array = np.round((u10-offset)/scale).astype('uint8')
    # band.append({'array':array, 'scale':scale, 'offset':offset,
    #              'description':'zonal wind speed', 'unittype':'m/s',
    #              'nodatavalue':255, 'parameter_range':[-25.4, 25.4]})
    # v10 = np.sin(winddirection)*windspeed
    # array = np.round((v10-offset)/scale).astype('uint8')
    # band.append({'array':array, 'scale':scale, 'offset':offset,
    #              'description':'meridional wind speed', 'unittype':'m/s',
    #              'nodatavalue':255, 'parameter_range':[-25.4, 25.4]})
    band = []
    indndv = np.where(land_mask == 1)
    windspeed = np.sqrt(u10**2 + v10**2)
    winddirection = np.mod(np.arctan2(v10, u10)*180./np.pi+360., 360.)
    offset, scale = vmin, (vmax-vmin)/254.
    np.clip(windspeed, vmin, vmax, out=windspeed)
    array = np.round((windspeed - offset) / scale).astype('uint8')
    array[indndv] = 255
    colortable = stfmt.format_colortable('noaa_wind', vmax=vmax, vmax_pal=vmax_pal,
                                         vmin=vmin, vmin_pal=vmin_pal)
    band.append({'array':array, 'scale':scale, 'offset':offset,
                 'description':'wind speed', 'unittype':'m/s',
                 'nodatavalue':255, 'parameter_range':[vmin, vmax],
                 'colortable':colortable})
    array = np.round(winddirection/360.*254.).astype('uint8')
    array[indndv] = 255
    band.append({'array':array, 'scale':360./254., 'offset':0.,
                 'description':'wind direction', 'unittype':'deg',
                 'nodatavalue':255, 'parameter_range':[0, 360.]})
    # Write geotiff
    if write_netcdf == False:
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)
    elif write_netcdf == True:
        print 'Write netcdf'
        # u/v -> bands
        band = []
        indndv = np.where(land_mask == 1)
        vmin = -vmax
        offset, scale = vmin, (vmax-vmin)/254.
        np.clip(u10, vmin, vmax, out=u10)
        array = np.round((u10 - offset) / scale).astype('uint8')
        array[indndv] = 255
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'wind u', 'unittype':'m s-1',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax],
                     'name':'u10m', 'long_name':'u component of horizontal wind',
                     'standard_name':'eastward_wind'})
        np.clip(v10, vmin, vmax, out=v10)
        array = np.round((v10 - offset) / scale).astype('uint8')
        array[indndv] = 255
        band.append({'array':array, 'scale':scale, 'offset':offset,
                     'description':'wind v', 'unittype':'m s-1',
                     'nodatavalue':255, 'parameter_range':[vmin, vmax],
                     'name':'v10m', 'long_name':'v component of horizontal wind',
                     'standard_name':'northward_wind'})
        # Write
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        metadata['spatial_resolution'] = 0.125 * 111000.
        dgcps = np.round(1. / np.array([0.125, 0.125])).astype('int')
        stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'grid_lonlat',
                           dgcps=dgcps)

