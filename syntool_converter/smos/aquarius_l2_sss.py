# -*- coding: utf-8 -*-

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

import netCDF4
import numpy
import os
import syntool_converter.utils.syntoolformat as stfmt
import datetime
import logging
from pyproj import Geod

logger = logging.getLogger(__name__)


def find_next_valid_lat(abs_lat, current, max_lat):
    """Return the index of the next absolute latitude which is below threshold.
    """
    l = len(abs_lat)
    i = current
    while (i < l) and (abs_lat[i] >= max_lat):
        i = i + 1
    return i


def compute_dlon(_lon, current):
    """Compute the difference between successive longitudes."""
    lon0 = _lon[current] - 180.0
    lon = numpy.mod(_lon - lon0, 360.0) + lon0
    dlon = lon[1:] - lon[:-1]
    return dlon


def split_by_valid_chunks(_lat, _lon, max_lat=85.0, offset=0):
    """Return swath slices that GDAL should be able to process without error.
    """
    abs_lat = numpy.abs(_lat)
    l = len(abs_lat)

    chunks = []
    i = find_next_valid_lat(abs_lat, 0, max_lat)
    if i >= l:
        return chunks
    chunk_start = i
    dlon = compute_dlon(_lon, chunk_start)
    last_dlon = None
    ended_with_invalid_lat = False
    while i < (l - 1):
        # Check latitude
        lat = abs_lat[i]
        if lat >= max_lat:
            logger.debug('Add chunk because lat above threshold: {}'.format(i))
            chunks.append(slice(offset + chunk_start, offset + i))

            i = find_next_valid_lat(abs_lat, i, max_lat)
            if i > (l - 1):
                ended_with_invalid_lat = True
                break
            chunk_start = i
            dlon = compute_dlon(_lon, chunk_start)
            last_dlon = None

        # Split to be able to increase GCPs density on high latitude only
        if 0 > (abs_lat[i] - 68.0) * (abs_lat[i+1] - 68.0):
            chunks.append(slice(offset + chunk_start, offset + i + 1))
            chunk_start = i  # i+1, but overlap added for smooth transition
            dlon = compute_dlon(_lon, chunk_start)
            last_dlon = None
            i = i + 1
            continue

        current_dlon = dlon[i]

        if 0 == current_dlon:
            # No longitudinal difference, proceed with next point
            i = i + 1
            continue

        # Check longitude
        if last_dlon is not None:
            lon_diff = numpy.mod(numpy.abs(current_dlon), 360.0)
            if  lon_diff > 90.0:
                logger.debug('Add chunk because lon gap too wide')
                chunks.append(slice(offset + chunk_start, offset + i + 1))
                chunk_start = i + 1
                dlon = compute_dlon(_lon, chunk_start)
                last_dlon = None
                i = i + 1
                continue
            elif numpy.sign(current_dlon) != last_dlon:
                logger.debug('Add chunk because lon not monotonous')
                chunks.append(slice(offset + chunk_start, offset + i + 1))
                chunk_start = i + 1
                dlon = compute_dlon(_lon, chunk_start)
                last_dlon = None
                i = i + 1
                continue

        last_dlon = numpy.sign(current_dlon)
        i = i + 1

    if not ended_with_invalid_lat:
        chunks.append(slice(offset + chunk_start,
                            offset + numpy.shape(_lat)[0]))
    return chunks


def process_slices(slices, dset, datagroup, beam, ref_dt, vmin, vmax, vmin_pal,
                   vmax_pal):
    """"""
    geod = Geod(ellps='WGS84')

    # Inner, middle and outer beam widths
    beam_widths = (94.0, 120.0, 156.0,)
    beam_width = beam_widths[beam]

    for k, chunk_rows in enumerate(slices):
        _chunk_lon = dset.variables['beam_clon'][chunk_rows, beam]
        chunk_size = numpy.shape(_chunk_lon)[0]

        if 0 >= chunk_size:
            # Empty chunk, skip
            continue

        chunk_lon0 = _chunk_lon[0] - 180.0
        chunk_lon = numpy.mod(_chunk_lon - chunk_lon0, 360.0) + chunk_lon0
        chunk_lat = dset.variables['beam_clat'][chunk_rows, beam]
        values = dset.variables['SSS'][chunk_rows, beam]
        values = numpy.ma.masked_where(values==-999., values)

        # Build GCPs
        dgcp = 32.
        if numpy.max(numpy.abs(chunk_lat)) > 75.0:
            dgcp = 4.

        ngcplin = numpy.ceil(chunk_lon.size / dgcp).astype('int32')
        _gcp_alongtrack = numpy.linspace(0, chunk_lon.size - 1, num=ngcplin)
        _gcp_indices = numpy.round(_gcp_alongtrack).astype('int32')
        _gcppix = numpy.array([-1.0, 0.5, 2.0])
        ngcppix = _gcppix.size
        gcppix = numpy.tile(_gcppix[numpy.newaxis, :],
                            (ngcplin, 1))
        gcpind = numpy.tile(_gcp_indices[:, numpy.newaxis],
                            (1, ngcppix)).astype('int32')
        gcplin = gcpind + 0.5

        # Compute swath direction
        _ind0 = numpy.minimum(gcpind, chunk_lon.size - 2)
        _ind1 = _ind0 + 1
        ind_same = numpy.where((chunk_lon[_ind0] == chunk_lon[_ind1]) &
                               (chunk_lat[_ind0] == chunk_lat[_ind1]))
        for ig_line, ig_pixel in zip(ind_same[0], ind_same[1]):
            if _ind1[ig_line, ig_pixel] < chunk_lon.size -1:
                _ind1[ig_line, ig_pixel] += 1
            else:
                _ind0[ig_line, ig_pixel] -= 1

        lat_diff = chunk_lat[_ind1] - chunk_lat[_ind0]
        lon_diff = chunk_lon[_ind1] - chunk_lon[_ind0]
        swath_dir = numpy.arctan2(lat_diff, lon_diff)

        # Compute GCPs geographical coordinates from the location of the beam
        # center and its width
        gcphei = numpy.zeros(gcppix.shape)
        gcplon = numpy.zeros(gcppix.shape)
        gcplat = numpy.zeros(gcppix.shape)
        for gcp_i in range(len(_gcp_indices)):
            ind = _gcp_indices[gcp_i]
            central_lon = chunk_lon[ind]
            central_lat = chunk_lat[ind]
            across_dir = -1 * numpy.rad2deg(numpy.mod(swath_dir[gcp_i][1],
                                            2 * numpy.pi))
            lon_a, lat_a, _ = geod.fwd(central_lon, central_lat,
                                       across_dir,
                                       1000.0 * 1.5 * beam_width)
            lon_b, lat_b, _ = geod.fwd(central_lon, central_lat,
                                       180.0 + across_dir,
                                       1000.0 * 1.5 * beam_width)

            gcplon[gcp_i][0] = lon_b
            gcplon[gcp_i][1] = central_lon
            gcplon[gcp_i][2] = lon_a

            gcplat[gcp_i][0] = lat_b
            gcplat[gcp_i][1] = central_lat
            gcplat[gcp_i][2] = lat_a

        # Fix longitudinal continuity
        half_ind = numpy.floor(len(_gcp_indices) * 0.5).astype('int32')
        gcplon0 = gcplon[half_ind, 1] - 180.0
        gcplon = numpy.mod(gcplon - gcplon0, 360.0) + gcplon0

        # Construct metadata/geolocation/band(s)
        sec = dset.variables['sec'][chunk_rows]
        start_dt = ref_dt + datetime.timedelta(seconds=sec[0])
        stop_dt = ref_dt + datetime.timedelta(seconds=sec[-1])
        dtime, time_range = stfmt.format_time_and_range(start_dt, stop_dt,
                                                        units='h')
        metadata = {}
        metadata['name'] = '{}_{}'.format(datagroup, k)
        metadata['time_range'] = time_range
        metadata['datetime'] = dtime
        metadata['datagroup'] = datagroup

        geolocation = {}
        geolocation['projection'] = stfmt.format_gdalprojection()
        geolocation['gcps'] = stfmt.format_gdalgcps(gcplon, gcplat, gcphei,
                                                    gcppix, gcplin)

        # Build mask
        land_frac = dset.variables['land_frac'][chunk_rows, beam]
        ice_frac = dset.variables['ice_frac'][chunk_rows, beam]
        scat_land_frac = dset.variables['scat_land_frac'][chunk_rows, beam]
        scat_ice_frac = dset.variables['scat_ice_frac'][chunk_rows, beam]
        mask = (values.mask |
                (land_frac > 0.1) |
                (ice_frac > 0.1))  # |
                # (scat_land_frac > 0.001) |
                # (scat_ice_frac > 0.001))

        # Pack data
        band = []
        offset, scale = vmin, (vmax - vmin) / 254.
        numpy.clip(values.data, vmin, vmax, out=values.data)
        array = numpy.round((values.data - offset) / scale).astype('uint8')
        array[numpy.where(mask)] = 255
        colortable = stfmt.format_colortable('matplotlib_jet',
                                             vmin=vmin, vmax=vmax,
                                             vmin_pal=vmin_pal,
                                             vmax_pal=vmax_pal)
        array = array[:, numpy.newaxis]
        band.append({'array':array,
                     'scale':scale,
                     'offset':offset,
                     'description':'sea surface salinity',
                     'unittype':'PSS',
                     'nodatavalue':255,
                     'parameter_range':[vmin, vmax],
                     'colortable':colortable})

        yield metadata, geolocation, band


def write_geotiff(metadata, geolocation, band, infile, outdir):
    """"""
    now = datetime.datetime.utcnow()

    metadata['product_name'] = 'AQUARIUS_L2_SSS'
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'Jet Propulsion Laboratory'
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(now)
    metadata['parameter'] = 'sea surface salinity'

    # Write geotiff
    tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
    stfmt.write_geotiff(tifffile, metadata, geolocation, band)


def aquarius_l2_sss(infile, outdir,
                vmin=31.825, vmax=38.175, vmin_pal=32, vmax_pal=38,
                write_netcdf=False):
    """"""
    # Read/Process data
    datagroup = os.path.splitext(os.path.basename(infile))[0]
    header, _ = datagroup.split('.', 1)
    _dt = datetime.datetime.strptime(header[1:], '%Y%j%H%M%S')
    _date = datetime.datetime(_dt.year, _dt.month, _dt.day)
    dset = netCDF4.Dataset(infile, 'r')
    time_start = dset['sec'][0]
    start_dt = _date + datetime.timedelta(seconds=time_start)
    time_stop = dset['sec'][-1]
    stop_dt = _date + datetime.timedelta(seconds=time_stop)
    _lat = dset['beam_clat'][:]
    _lon = dset['beam_clon'][:]

    beam_names = ('inner', 'middle', 'outer',)

    # Process ascending and descending phases separately to avoid issues with
    # GDAL reprojection when the longitudinal extent of the shape is too big.
    half_split = numpy.ceil(len(dset['sec']) / 2.0).astype('int32')
    for beam in range(3):
        beam_lon = _lon[:half_split, beam]
        beam_lat = _lat[:half_split, beam]
        slices = split_by_valid_chunks(beam_lat, beam_lon, 85.4, 0)
        _datagroup = '{}_p0_{}'.format(datagroup, beam_names[beam])
        for meta, geoloc, data in process_slices(slices, dset, _datagroup,
                                                 beam, _date, vmin, vmax,
                                                 vmin_pal, vmax_pal):
            write_geotiff(meta, geoloc, data, infile, outdir)

        beam_lon = _lon[half_split:, beam]
        beam_lat = _lat[half_split:, beam]
        slices = split_by_valid_chunks(beam_lat, beam_lon, 85.4, half_split)
        _datagroup = '{}_p1_{}'.format(datagroup, beam_names[beam])
        for meta, geoloc, data in process_slices(slices, dset, _datagroup,
                                                 beam, _date, vmin, vmax,
                                                 vmin_pal, vmax_pal):
            write_geotiff(meta, geoloc, data, infile, outdir)
    dset.close()
