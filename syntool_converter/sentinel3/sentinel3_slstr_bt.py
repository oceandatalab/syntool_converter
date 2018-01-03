# -*- encoding: utf-8 -*-

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

import os
import sys
import numpy
import logging
import syntool_converter.utils.syntoolformat as stfmt
from scipy.misc import bytescale
from datetime import datetime
import syntool_converter.sentinel3.slstr as slstr

logger = logging.getLogger(__name__)

# https://earth.esa.int/web/sentinel/technical-guides/sentinel-3-slstr/level-1/vis-and-swir-radiances
# https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-slstr/data-formats/level-1
# https://earth.esa.int/web/sentinel/user-guides/sentinel-3-slstr/product-types/level-1b
# https://sentinel.esa.int/documents/247904/685236/Sentinel-3_User_Handbook#%5B%7B%22num%22%3A235%2C%22gen%22%3A0%7D%2C%7B%22name%22%3A%22XYZ%22%7D%2C54%2C628%2C0%5D
# S1 555   nm            B
# S2 659   nm            G
# S3 865   nm            R
# S4 1.375 microm
# S5 1.610 microm
# S6 2.25  microm
# S7 3.74  microm        IR
# S8 10.85 microm        IR
# S9 12    microm        IR


default_minmax = {'S7': (270, 305),
                  'S8': (270, 305),
                  'S9': (270, 305)}


def build_mask_ir(channels, quality_flags, raw_cloud_flags, tie_lat, lat_crop):
    """Build masks to compute contrast and to mask invalid data"""

    bits = {'coastline': 0,
            'ocean': 1,
            'tidal': 2,
            'land': 3,
            'inland_water': 4,
            'unfilled': 5,
            'spare': 6,
            'spare_': 7,
            'cosmetic': 8,
            'duplicate': 9,
            'day': 10,
            'twilight': 11,
            'sun_glint': 12,
            'snow': 13,
            'summary_cloud': 14,
            'summary_pointing': 15}

    cloud_bits = {'visible': 0,
                  'threshold': 1,
                  'small_histogram1.6': 2,
                  'large_histogram1.6': 3,
                  'small_histogram2.25': 4,
                  'large_histogram2.25': 5,
                  'spatial_coherence': 6,
                  'gross_cloud': 7,
                  'thin_cirrus': 8,
                  'medium_high': 9,
                  'fog_low_stratus': 10,
                  '11_12_view_difference': 11,
                  '3.7_11_view_difference': 12,
                  'thermal_histogram': 13}

    # Invalidate if these flags are all off
    on_flags = numpy.ushort(0)
    on_flags = on_flags + numpy.ushort(1 << bits['day'])

    # Invalidate if any of these flags is on
    off_flags = numpy.ushort(0)
    off_flags = off_flags + numpy.ushort(1 << bits['coastline'])
    off_flags = off_flags + numpy.ushort(1 << bits['land'])
    off_flags = off_flags + numpy.ushort(1 << bits['snow'])

    # Intialize mask to compute histograms on selected valid values
    contrast_mask = numpy.zeros(quality_flags.shape, dtype='bool')

    # Initialize mask to remove invalid data
    data_mask = numpy.zeros(quality_flags.shape, dtype='bool')

    lat_mask = (numpy.abs(tie_lat) > lat_crop)
    contrast_mask = (contrast_mask | lat_mask)
    data_mask = (data_mask | lat_mask)

    land_mask = slstr.build_land_mask(quality_flags, bits)

    cloud_mask = slstr.build_cloud_mask(raw_cloud_flags, cloud_bits)

    # Do not mask data at night or twilight for ir channels.
    # Apply land mask only for contrast computation: land data must still
    # be visible in the output.
    contrast_mask = (contrast_mask |
                     cloud_mask |
                     land_mask |
                     (numpy.bitwise_and(quality_flags, off_flags) > 0))

    data_mask = (data_mask |
                 cloud_mask |
                 (numpy.bitwise_and(quality_flags, off_flags) > 0))

    return contrast_mask, data_mask


def get_valid_data_ir(tie_lon, tie_lat, file_range, band,
                      mask, bandname, month):
    """ Build valid data mask from to detect range value."""

    extra_data_mask = None
    updated_min = None
    # Compute min from climatology
    list_min = []
    list_max = []

    if file_range is not None:
        # Find the range file which corresponds to the time coverage of the
        # granule
        base, ext = os.path.splitext(file_range)
        static_base, _ = base.rsplit('_', 1)
        range_file = '{}_{:02d}{}'.format(static_base, month, ext)
        if not os.path.exists(range_file):
            logger.error('{} not found'.format(range_file))
            sys.exit(1)

        with open(range_file, 'r') as fp:
            lines = fp.readlines()

        if 64800 != len(lines):  # 64800 = 360 * 180
            logger.error('The range file should have a 1 degree resolution '
                         '(64800 lines)')
            sys.exit(1)

        # Find longitude range
        lons180 = numpy.mod(tie_lon + 180.0, 360.0) - 180.0
        lon_min_int = numpy.floor(numpy.min(lons180)).astype('int32')
        lon_max_int = numpy.floor(numpy.max(lons180)).astype('int32')
        range_lon = range(lon_min_int, lon_max_int + 1)  # [lon_min, lon_max]

        # Use mean distance from the Greenwich Meridian to determine if the
        # granule overlaps the International Date Line.
        over_idl = bool(((0 > lon_min_int * lon_max_int) and
                    (90.0 < numpy.mean(numpy.abs(lons180)))))
        if over_idl is True:
            # Build the longitudes range starting with positive values, then
            # then negative ones to fix continuity issues
            lon_min = numpy.max(lons180[lons180 < 0.0])
            lon_max = numpy.min(lons180[lons180 > 0.0])
            lon_min_int = numpy.floor(lon_min).astype('int32')
            lon_max_int = numpy.floor(lon_max).astype('int32')
            range_lon = range(lon_max_int, 180)  # [lon_max, 180[
            range_lon.extend(range(-180, lon_min_int + 1))  # [-180, lon_min]

        # Find latitude range
        lat_min_int = numpy.floor(numpy.min(tie_lat)).astype('int32')
        lat_max_int = numpy.floor(numpy.max(tie_lat)).astype('int32')
        range_lat= range(lat_min_int, lat_max_int + 1)  # [lat_min, lat_max]

        # Extract min/max for each value in the (range_lon x range_lat) domain
        # The "ll" prefix stands for "lower left"
        for _lllat in range_lat:
            lllat = _lllat + 90  # express latitude in [0, 180]
            for _lllon in range_lon:
                lllon = (_lllon + 360) % 360  # express longitude in [0, 360[

                line_ind = lllat * 360 + lllon
                iline = lines[line_ind].split(',')
                tmpmin = float(iline[4])
                if not numpy.isnan(tmpmin):
                    list_min.append(tmpmin)
                tmpmax = float(iline[5])
                if not numpy.isnan(tmpmax):
                    list_max.append(tmpmax)

    vmin_mask, vmax_mask = default_minmax[bandname]
    if len(list_min) > 0:
        vmin_mask = min(list_min)
        # vmax_mask = max(list_max)

    # Constrain min value for upwelling areas and use a threshold for salted
    # water (earth TPNC). SUIM
    lat_correction = 3
    value_min = max(vmin_mask - lat_correction, 269.35)

    # A value which is below the min obtained with the climatology
    # probably belongs to a cloud, mask it
    min_mask = numpy.ma.masked_less(band, value_min)

    valid_data = band.data[~mask & ~min_mask.mask]
    extra_data_mask = min_mask.mask
    updated_min = value_min

    return valid_data, extra_data_mask, updated_min


def sentinel3_slstr_bt(infile, outdir, vmin=None, vmax=None,
                       min_percentile=2.0, channels='ir', file_range=None,
                       write_netcdf=False, log_path=None, lat_crop=80.0):

    t0 = datetime.utcnow()
    # Process nadir data
    view = 'n'
    fname = 'BT'
    sltype = 'i'
    if type(channels) is list or type(channels) is tuple:
        bandnames = channels
        product_name = 'Sentinel-3_SLSTR'
    elif 'ir' == channels:
        bandnames = ('S8',)
        product_name = 'Sentinel-3_SLSTR_IR'
    else:
        raise Exception('channels must be either "ir", or a tuple of band')

    # Read coordinates and compute gcps
    (syntool_stats, metadata, geolocation, tie_lon, tie_lat, slice_lat0,
        slice_lat1, __, __, ngcps, month) = slstr.read_geometry(infile,
                                             bandnames, fname, sltype, view,
                                             product_name, vmin, vmax,
                                             log_path)

    # Compute masks
    logger.info('Build masks')
    t_start = datetime.utcnow()
    quality_flags = slstr.read_mask(infile, sltype, view, slice_lat0,
                                    slice_lat1)
    raw_cloud_flags = slstr.read_cloud_mask(infile, sltype, view, slice_lat0,
                                            slice_lat1)
    contrast_mask, data_mask = build_mask_ir(channels, quality_flags,
                                             raw_cloud_flags, tie_lat,
                                             lat_crop)
    t_stop = datetime.utcnow()
    syntool_stats['mask_computation'] = (t_stop - t_start).total_seconds()

    # Read band to compute histograms
    logger.info('Construct bands')
    t_start = datetime.utcnow()

    bands = []
    # Initialize min and max values
    if vmin is None:
        vmin = [None] * len(bandnames)
    if vmax is None:
        vmax = [None] * len(bandnames)
    _vmin = list(vmin)
    _vmax = list(vmax)
    for band_index in range(len(bandnames)):
        bandname = bandnames[band_index]
        fieldname = slstr.get_field_name(fname, bandname, sltype, view)
        band = slstr.read_band(infile, bandname, fieldname, slice_lat0,
                               slice_lat1)

        logger.info('\tSet contrast')
        valid_ratio_lower_threshold = 0.001  # 0.1%

        # Select valid data to compute histograms
        valid_data_mask = (band.mask | contrast_mask)
        valid_data, extra_data_mask, updated_min = get_valid_data_ir(tie_lon,
                                                   tie_lat, file_range, band,
                                                   valid_data_mask, bandname,
                                                   month)

        if extra_data_mask is not None:
            data_mask = (data_mask | extra_data_mask)
        if updated_min is not None:
            _vmin = [updated_min, ]
            _min = updated_min

        # No need to produce an output if all data values are masked
        if numpy.all(data_mask):
            logger.warn('No valid value found for band {}'.format(bandname))
            sys.exit(0)

        # Retrieve minimum and maximum values from default or valid_data
        # histograms
        valid_ratio = float(valid_data.size) / float(band.data.size)
        syntool_stats[bandname]['valid_ratio'] = valid_ratio
        if valid_ratio_lower_threshold >= valid_ratio:
            _min, _max = slstr.apply_default_min_max(default_minmax, bandname,
                                                     _vmin[band_index],
                                                     _vmax[band_index],
                                                     syntool_stats)

        else:
            _min, _max = slstr.fromband_min_max(valid_data, bandname,
                                                _vmin[band_index],
                                                _vmax[band_index],
                                                syntool_stats,
                                                min_percentile=min_percentile,
                                                max_percentile=99.99)

        _vmin[band_index] = _min
        _vmax[band_index] = _max
        logger.info('\tContrast : vmin={} / vmax={}'.format(_vmin[band_index],
                                                            _vmax[band_index]))
    min_values = [_vmin[band_index] for band_index in range(len(bandnames))]
    max_values = [_vmax[band_index] for band_index in range(len(bandnames))]

    t_stop = datetime.utcnow()
    syntool_stats['minmax_computation'] = (t_stop - t_start).total_seconds()
    syntool_stats['final_min'] = float(numpy.min(min_values))
    syntool_stats['final_max'] = float(numpy.max(max_values))

    _min = numpy.min(min_values)
    _max = numpy.max(max_values)
    scale = (_max - _min) / 254.
    offset = _min
    # Construct bands
    for band_index in range(len(bandnames)):
        bandname = bandnames[band_index]
        fieldname = slstr.get_field_name(fname, bandname, sltype, view)
        band = slstr.read_band(infile, bandname, fieldname, slice_lat0,
                               slice_lat1)

        bnd = band.data

        logger.info('\tBytescaling')
        byte = bytescale(bnd, cmin=_min, cmax=_max, low=0, high=254)
        description = '{} {} (log)'.format(bandname, fname)
        if band.mask is not numpy.ma.nomask:
            byte[band.mask] = 255

        # mask night data for rgb and invalid data for ir (cloud, land,
        # range value). Also mask data for extreme latitudes
        byte[data_mask] = 255

        band_range = [_vmin[band_index], _vmax[band_index]]
        description = '{} {}'.format(bandname, fname)  # no log for IR
        colortable = stfmt.format_colortable('cerbere_medspiration',
                                             vmax=_max, vmax_pal=_max,
                                             vmin=_min, vmin_pal=_min)
        bands.append({'array': byte,
                      'plot': band.data,
                      'scale': scale,
                      'offset': offset,
                      'description': description,
                      'unittype': '',
                      'nodatavalue': 255,
                      'parameter_range': band_range,
                      'colortable': colortable})

        if write_netcdf:
            bands[-1]['name'] = bandname
            bands[-1]['long_name'] = bandname
            bands[-1]['unittype'] = '1'

    logger.info('Make sure nodata are at the same locations in all bands')
    mask = numpy.any([_band['array'] == 255 for _band in bands], axis=0)
    for band in bands:
        band['array'][mask] = 255

    t_stop = datetime.utcnow()
    syntool_stats['bytescaling'] = (t_stop - t_start).total_seconds()

    if write_netcdf:
        metadata['spatial_resolution'] = 1000
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        stfmt.write_netcdf(ncfile, metadata, geolocation, bands, 'swath',
                           ngcps=ngcps)
    else:
        logger.info('Write geotiff')
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, bands)

    logger.info(datetime.utcnow() - t0)
    syntool_stats['total_time'] = (datetime.utcnow() - t0).total_seconds()
    if log_path is not None:
        import json
        full_path = os.path.normpath(infile)
        file_path = os.path.basename(full_path)
        file_name, _ = os.path.splitext(file_path)
        stats_path = os.path.join(log_path, '{}.json'.format(file_name))
        with open(stats_path, 'w') as f:
            json.dump(syntool_stats, f)
