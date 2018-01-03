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

from __future__ import print_function
from netCDF4 import Dataset
import numpy as np
import numpy.ma
from datetime import datetime
import syntool_converter.utils.syntoolformat as stfmt
import logging
import os

logger = logging.getLogger(__name__)

def make_continuous_lon(lon, max_gap):
    """Returns continuous longitudes and indices for rows where the difference
       with the previous row equals or exceeds max_gap degrees."""
    irow = 0
    nrow = lon.shape[0]

    # Enforce continuity along track (the longitudinal extent may cover more
    # than 360 degrees)
    dlon = lon[1:, :] - lon[:-1,:]
    dlon[dlon > 180.] = dlon[dlon > 180.] - 360.
    dlon[dlon < -180.] = dlon[dlon < -180.] + 360.
    for irow in xrange(1, nrow):
        prev_row = irow - 1
        lon[irow, :] = lon[prev_row, :] + dlon[prev_row, :]
    extra_large_delta = np.where(max_gap <= np.abs(dlon))

    # Enforce continuity across-track (more than 360 degrees is not possible)
    for irow in xrange(0, nrow):
        lon0 = lon[irow,0] + 180.
        lon[irow,:] = np.mod(lon[irow,:] - lon0, 360.) + lon0

    # Add a +1 offset since there is no dlon applies to lon[1:]
    rows_with_extra_large_delta = 1 + np.unique(extra_large_delta[0])
    return lon, rows_with_extra_large_delta.tolist()


def split_lon_chunks(lon, max_size=180.):
    """"""
    nrows = lon.shape[0]
    slices = []
    slice_start = 0
    slice_stop = 0
    ref = np.max(lon[slice_start, :])
    logger.debug('Split longitudes')
    logger.debug('\tref: {}'.format(ref))

    for i in range(0, nrows):
        if max_size > np.abs(ref - np.min(lon[i, :])):
            slice_stop += 1
        else:
            slices.append(slice(slice_start, slice_stop))
            slice_start = i
            slice_stop = i
            ref = np.max(lon[slice_start, :])
            logger.debug('\tref: {}'.format(ref))
    if slice_stop > slices[-1].stop:
        slices.append(slice(slice_start, slice_stop))
    return slices

def _to_geotiff(infile, outdir, vmin, vmax, vmin_pal, vmax_pal,
                nrow, ncell, dtime, time_range, datagroup, platform, dset):
    """ """
    # It is necessary to create a 1-row overlap between the two half-orbits in
    # order for the merger to build a continuous shape
    extra_row = 1

    side_slices = {'left': slice(0, ncell / 2),
                   'right': slice(ncell / 2, ncell)}

    for side in side_slices:
        side_slice = side_slices[side]

        side_lat = dset.variables['lat'][:, side_slice]
        side_lon = dset.variables['lon'][:, side_slice]

        # Make longitudes continuous but detect gaps bigger than 90 degrres:
        # it will be necessary to split the dataset where these gaps appear
        # because gdal may not be able to interpolate GCPs correctly later
        # in the ingestor (it says the shape intersects itself).
        side_lon, splits = make_continuous_lon(side_lon, max_gap=60.)

        # Add begin and end to the splits list
        splits.append(0)
        splits.append(nrow)

        # Add an artificial split in the middle of the swath
        half_split = nrow/2
        splits.append(half_split)

        # Make sure there are no duplicates and sort the indices in ascending
        # order
        splits = list(set(splits))
        splits.sort()

        l = len(splits)
        parts = []
        extra_row_part = -1
        for i in range(l-1):
            extra_row = 0
            if splits[i+1] == half_split:
                # Since the half split does not correspond to a void area in
                # the original data, add an extra row so that there will be no
                # visible cut when the ingestor merges the parts back together
                extra_row_part = i
                extra_row = 1
            part = slice(splits[i], extra_row + splits[i+1])
            parts.append(part)
        nparts = len(parts)

        last_lon = None
        last_lat = None
        for part_no in range(0, nparts):
            dataset_name = '{}_{}_{}'.format(datagroup, side, part_no)

            part_slice = parts[part_no]
            swath_slice = [part_slice, side_slice]

            lat = side_lat[part_slice]
            lon = side_lon[part_slice]

            # Make sure there are no longitudes < -180 or the rastertiles
            # plugin of syntool-ingestor will not produce anything
            while 180. < np.min(lon):
                lon = lon - 360.
            while -180. > np.min(lon):
                lon = lon + 360.

            # Extract wind speed (module)
            wind_speed = dset.variables['wind_speed'][swath_slice]

            # Extract wind direction and make it counter-clockwise from east
            # (as expected by the ingestor)
            wind_dir = dset.variables['wind_dir'][swath_slice]
            wind_dir = np.mod(90. - wind_dir, 360.)

            # Build GCPs
            dgcp = 16.
            ngcps = (np.ceil(np.array(lon.shape) / dgcp) + 1.).astype('int32')
            pix = np.linspace(0, lon.shape[1] - 1, num=ngcps[1]).round()
            lin = np.linspace(0, lon.shape[0] - 1, num=ngcps[0]).round()
            pix2d, lin2d = np.meshgrid(pix.astype('int32'),
                                       lin.astype('int32'))
            gcplon = lon[lin2d, pix2d]
            gcplat = lat[lin2d, pix2d]
            gcppix = pix2d + 0.5
            gcplin = lin2d + 0.5
            gcphei = np.zeros(ngcps)

            if part_no == extra_row_part + 1:
                # Overwrite first GCP, but first you must make sure that
                # longitudes are in the same -180/+180 range.
                lon0 = gcplon[1, 0] - 180.
                last_lon[:] = np.mod(last_lon[:] - lon0, 360.) + lon0
                gcplon[0] = last_lon
                gcplat[0] = last_lat

            if part_no == extra_row_part:
                # Mask last row since it also exists in the next part
                wind_speed[-extra_row, :] = numpy.ma.masked
                wind_dir[-extra_row, :] = numpy.ma.masked
                # Store last GCPs lat/lon
                extra_row = 1
                last_lon = gcplon[-extra_row]
                last_lat = gcplat[-extra_row]

            # Construct metadata/geolocation/band(s)
            print('Construct metadata/geolocation/band(s)')
            utcnow = datetime.utcnow()
            metadata = {}
            metadata['product_name'] = '{}_ASCAT_L2B'.format(platform)
            metadata['datagroup'] = datagroup
            metadata['name'] = dataset_name
            metadata['datetime'] = dtime
            metadata['time_range'] = time_range
            metadata['source_URI'] = infile
            metadata['conversion_software'] = 'Syntool'
            metadata['conversion_version'] = '0.0.0'
            metadata['conversion_datetime'] = stfmt.format_time(utcnow)
            metadata['parameter'] = ['wind speed', 'wind direction']
            metadata['original_x_size'] = ncell
            metadata['original_y_size'] = nrow
            geolocation = {}
            geolocation['projection'] = stfmt.format_gdalprojection()
            geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                        gcppix, gcplin)
            print('Write geotiff')

            band = []
            offset, scale = vmin, (vmax-vmin)/254.
            clipped = np.clip(np.ma.getdata(wind_speed), vmin, vmax)
            array = np.round((clipped - offset) / scale).astype('uint8')
            array[np.ma.getmaskarray(wind_speed)] = 255
            colortable = stfmt.format_colortable('noaa_wind',
                                                 vmax=vmax, vmax_pal=vmax_pal,
                                                 vmin=vmin, vmin_pal=vmin_pal)
            band.append({'array':array, 'scale':scale, 'offset':offset,
                         'description':'wind speed', 'unittype':'m/s',
                         'nodatavalue':255, 'parameter_range':[vmin, vmax],
                         'colortable':colortable})
            clipped = np.clip(np.ma.getdata(wind_dir), 0, 360)
            array = np.round(clipped / 360. * 254.).astype('uint8')
            array[np.ma.getmaskarray(wind_dir)] = 255
            band.append({'array':array, 'scale':360./254., 'offset':0.,
                         'description':'wind direction', 'unittype':'deg',
                         'nodatavalue':255, 'parameter_range':[0, 360.]})

            tifffile = stfmt.format_tifffilename(outdir, metadata,
                                                 create_dir=True)
            stfmt.write_geotiff(tifffile, metadata, geolocation, band)


def ascat_l2b(infile, outdir,
              vmin=0., vmax=25.4, vmin_pal=0., vmax_pal=50*0.514,
              write_netcdf=False):
    """
    """
    # Read/Process data
    dset = Dataset(infile)
    nrow = len(dset.dimensions['NUMROWS'])
    ncell = len(dset.dimensions['NUMCELLS'])
    if ncell != 82:
        raise Exception('Expects NUMCELLS=82 (KNMI ASCAT L2B 12.5km).')
    source = dset.source
    if 'metop-a' in source.lower():
        platform = 'Metop-A'
    elif 'metop-b' in source.lower():
        platform = 'Metop-B'
    else:
        raise Exception('Platform ?')
    start_time = datetime.strptime(dset.start_date + dset.start_time,
                                   '%Y-%m-%d%H:%M:%S')
    stop_time = datetime.strptime(dset.stop_date + dset.stop_time,
                                  '%Y-%m-%d%H:%M:%S')
    (dtime, time_range) = stfmt.format_time_and_range(start_time, stop_time,
                                                      units='ms')
    datagroup = os.path.splitext(os.path.basename(infile))[0]
    if write_netcdf == False:
        # Create GeoTIFF file
        _to_geotiff(infile, outdir, vmin, vmax, vmin_pal, vmax_pal,
                    nrow, ncell, dtime, time_range, datagroup, platform, dset)
    else:
        for i in range(2):
            if i == 0: # left swath
                swath_slice = [slice(0, nrow), slice(0, ncell / 2)]
                dataset_name = '{}_left'.format(datagroup)
            else: # right swath
                swath_slice = [slice(0, nrow), slice(ncell / 2, ncell)]
                dataset_name = '{}_right'.format(datagroup)
            lat = dset.variables['lat'][swath_slice]
            lon = dset.variables['lon'][swath_slice]
            irow = 0
            while irow < nrow:
                lon0 = lon[irow, 0] - 180
                lon[irow:, :] = np.mod(lon[irow:, :] - lon0, 360) + lon0
                notcont = (lon[irow + 1:, :] < lon0 + 90) & (lon[irow:-1, :] > lon0 + 270) | \
                          (lon[irow + 1:, :] > lon0 + 270) & (lon[irow:-1, :] < lon0 + 90)
                indnotcont = np.where(notcont.any(axis=1))[0]
                if indnotcont.size == 0:
                    irow = nrow
                else:
                    indnotcont = indnotcont.min()
                    if indnotcont == 0:
                        raise Exception('Unexpected longitudes.')
                    irow = irow + indnotcont
            ind = np.where(np.abs(lon[1:, :] - lon[:-1, :]) > 180.)
            if ind[0].size != 0:
                raise Exception('Failed to make longitudes continuous.')
            if lon[nrow / 2, ncell / 4] > 180:
                lon -= 360
            elif lon[nrow / 2, ncell / 4] < -180:
                lon += 360

            print(lon)

            wind_speed = dset.variables['wind_speed'][swath_slice]
            wind_dir = dset.variables['wind_dir'][swath_slice]
            wind_dir = np.mod(90. - wind_dir, 360.)
            dgcp = 16.
            ngcps = (np.ceil(np.array(lon.shape) / dgcp) + 1.).astype('int32')
            pix = np.linspace(0, lon.shape[1] - 1, num=ngcps[1]).round().astype('int32')
            lin = np.linspace(0, lon.shape[0] - 1, num=ngcps[0]).round().astype('int32')
            pix2d, lin2d = np.meshgrid(pix, lin)
            gcplon = lon[lin2d, pix2d]
            gcplat = lat[lin2d, pix2d]
            gcppix = pix2d + 0.5
            gcplin = lin2d + 0.5
            gcphei = np.zeros(ngcps)
            # Construct metadata/geolocation/band(s)
            print('Construct metadata/geolocation/band(s)')
            metadata = {}
            metadata['product_name'] = '{}_ASCAT_L2B'.format(platform)
            metadata['datagroup'] = datagroup
            metadata['name'] = dataset_name
            metadata['datetime'] = dtime
            metadata['time_range'] = time_range
            metadata['source_URI'] = infile
            metadata['conversion_software'] = 'Syntool'
            metadata['conversion_version'] = '0.0.0'
            metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
            metadata['parameter'] = ['wind speed', 'wind direction']
            geolocation = {}
            geolocation['projection'] = stfmt.format_gdalprojection()
            geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                        gcppix, gcplin)
            print('Write netcdf')
            # u/v -> bands
            band = []
            u = wind_speed * np.cos(np.deg2rad(wind_dir))
            v = wind_speed * np.sin(np.deg2rad(wind_dir))
            mask = np.ma.getmaskarray(u) | np.ma.getmaskarray(v)
            vmin = -vmax
            offset, scale = vmin, (vmax-vmin)/254.
            clipped = np.clip(np.ma.getdata(u), vmin, vmax)
            array = np.round((clipped - offset) / scale).astype('uint8')
            array[mask] = 255
            band.append({'array':array, 'scale':scale, 'offset':offset,
                         'description':'wind u', 'unittype':'m s-1',
                         'nodatavalue':255, 'parameter_range':[vmin, vmax],
                         'name':'u', 'standard_name':'eastward_wind'})
            clipped = np.clip(np.ma.getdata(v), vmin, vmax)
            array = np.round((clipped - offset) / scale).astype('uint8')
            array[mask] = 255
            band.append({'array':array, 'scale':scale, 'offset':offset,
                         'description':'wind v', 'unittype':'m s-1',
                         'nodatavalue':255, 'parameter_range':[vmin, vmax],
                         'name':'v', 'standard_name':'northward_wind'})
            # Write
            ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
            metadata['spatial_resolution'] = 12500.
            stfmt.write_netcdf(ncfile, metadata, geolocation, band, 'swath',
                               ngcps=gcplon.shape)

