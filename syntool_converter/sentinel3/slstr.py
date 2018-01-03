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
import errno
import numpy
import netCDF4
import logging
import syntool_converter.utils.syntoolformat as stfmt
from datetime import datetime
from scipy.ndimage.morphology import binary_dilation

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


def build_land_mask(quality_flags, bits):
    """ Extand land mask to avoid detecting extreme values near the coast
    while detecting min and max values. """
    land_mask = numpy.zeros(quality_flags.shape, dtype='bool')
    land_bit = numpy.ushort(1 << bits['land'])
    coast_bit = numpy.ushort(1 << bits['coastline'])
    land_mask = (land_mask |
                 (numpy.bitwise_and(quality_flags, land_bit) > 0) |
                 (numpy.bitwise_and(quality_flags, coast_bit) > 0))

    # Dilatation of land mask to remove coastal values
    dil_kern = numpy.ones((5, 5), dtype='bool')
    land_mask = binary_dilation(land_mask, structure=dil_kern)

    return land_mask


def build_cloud_mask(raw_cloud_flags, cloud_bits):
    """ Construct a cloud mask from the raw flags in the SLSTR flag mask. Bits
    to consider in the mask have been selected to filter as much cloud as
    possible without removing frontal structures."""
    bits_to_mask = ('visible',
                    'threshold',
                    'small_histogram1.6',
                    # 'large_histogram1.6',
                    'small_histogram2.25',
                    'gross_cloud',
                    'thin_cirrus',
                    'medium_high',
                    'fog_low_stratus',
                    # '3.7_11_view_difference',
                    '11_12_view_difference')
    mask_ref = numpy.ushort(0)
    # Enable selected mask bits
    for mask_name in bits_to_mask:
        mask_ref = mask_ref + (1 << cloud_bits[mask_name])

    # Apply spatial coherence only if thermal histogram is flagged to
    # avoid masking fronts
    # 1 - Compute thermal_histogram mask
    unmask_ref = numpy.ushort(0)
    unmask_ref = unmask_ref + (1 << cloud_bits['thermal_histogram'])
    # 2 - Compute spatial coherence mask
    spatial = numpy.ushort(0)
    spatial = spatial + (1 << cloud_bits['spatial_coherence'])
    # 3 - Combine masks
    spatial_mask = ((numpy.bitwise_and(raw_cloud_flags, unmask_ref) > 0) &
                    (numpy.bitwise_and(raw_cloud_flags, spatial) > 0))

    # Build cloud mask
    cloud_mask = numpy.zeros(shape=raw_cloud_flags.shape, dtype='bool')
    cloud_mask = (cloud_mask |
                  spatial_mask |
                  (numpy.bitwise_and(raw_cloud_flags, mask_ref) > 1))

    return cloud_mask


def get_field_name(fname, bandname, sltype, view):
    """ Construct fieldname from name of variable, band, type of sensor and
    view. """
    fieldname = '{}_{}_{}{}'.format(bandname, fname, sltype, view)
    return fieldname


def read_geometry(infile, bandnames, fname, sltype, view, product_name, vmin,
                  vmax, log_path):
    """ Read coordinates from geometry SLSTR file """
    # Initiate log file
    syntool_stats = {}
    for bandname in bandnames:
        syntool_stats[bandname] = {}
    if log_path is not None and not os.path.exists(log_path):
        try:
            os.makedirs(log_path)
        except OSError:
            _, e, _ = sys.exc_info()
            if e.errno != errno.EEXIST:
                raise

    if vmin is None:
        vmin = [None] * len(bandnames)
    if vmax is None:
        vmax = [None] * len(bandnames)

    full_path = os.path.normpath(infile)
    file_path = os.path.basename(full_path)
    file_name, _ = os.path.splitext(file_path)

    geo_filename = 'geodetic_{}{}.nc'.format(sltype, view)
    geo_path = os.path.join(infile, geo_filename)
    lat_varname = 'latitude_{}{}'.format(sltype, view)
    lon_varname = 'longitude_{}{}'.format(sltype, view)

    # Extract geo coordinates information
    geo_handler = netCDF4.Dataset(geo_path, 'r')
    nrow = geo_handler.dimensions['rows'].size
    nrow_all = nrow
    ncell = geo_handler.dimensions['columns'].size
    ncell_all = ncell
    lon = geo_handler.variables[lon_varname][:]
    tie_lon = numpy.ma.array(lon)
    tie_lon._sharedmask = False
    lat = geo_handler.variables[lat_varname][:]
    tie_lat = numpy.ma.array(lat)
    start_time_str = geo_handler.start_time
    stop_time_str = geo_handler.stop_time
    geo_handler.close()

    # Handle longitude range
    dlon = tie_lon[1:, :] - tie_lon[:-1, :]
    if 180.0 <= numpy.max(numpy.abs(dlon)):
        tie_lon[tie_lon < 0.0] = tie_lon[tie_lon < 0.0] + 360.0

    # Extract time coordinates information
    # Format time information
    start_time = datetime.strptime(start_time_str,
                                   "%Y-%m-%dT%H:%M:%S.%fZ")
    end_time = datetime.strptime(stop_time_str,
                                 "%Y-%m-%dT%H:%M:%S.%fZ")
    (dtime, time_range) = stfmt.format_time_and_range(start_time, end_time,
                                                      units='ms')
    month = start_time.month

    parameters = ['{} TOA {}'.format(bnd, fname) for bnd in bandnames]
    metadata = {}
    metadata['product_name'] = '{}'.format(product_name)
    metadata['name'] = file_name
    metadata['datetime'] = dtime
    metadata['time_range'] = time_range
    metadata['source_URI'] = infile
    metadata['source_provider'] = 'ESA'
    metadata['processing_center'] = 'OceanDataLab'
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    metadata['parameter'] = parameters
    metadata['type'] = 'remote sensing'
    metadata['sensor_type'] = 'dual scan temperature radiometer'
    metadata['sensor_name'] = 'SLSTR'
    metadata['sensor_platform'] = 'Sentinel-3'

    # Crop high latitude to avoid projection issues
    LAT_MAX = 89.0
    ind_valid_cols = numpy.where(numpy.abs(tie_lat).max(axis=0) <= LAT_MAX)[0]
    slice_lat0 = slice(None)
    slice_lat1 = slice(numpy.min(ind_valid_cols),
                       numpy.max(ind_valid_cols) + 1)
    tie_lat = tie_lat[slice_lat0, slice_lat1]
    tie_lon = tie_lon[slice_lat0, slice_lat1]
    nrow, ncell = tie_lat.shape

    # Handle longitude continuity
    dlon = tie_lon[1:, :] - tie_lon[:-1, :]
    if 180.0 <= numpy.max(numpy.abs(dlon)):
        lon0 = tie_lon[0, 0] + 180.0
        tie_lon._sharedmask = False
        tie_lon[:, :] = numpy.mod(tie_lon[:, :] - lon0, 360.0) + lon0

    # Compute GCPs
    tie_row = numpy.linspace(0, nrow - 1, num=tie_lon.shape[0])
    tie_cell = numpy.linspace(0, ncell - 1, num=tie_lon.shape[1])
    tie_facrow = (nrow - 1.) / (tie_lon.shape[0] - 1.)
    tie_faccell = (ncell - 1.) / (tie_lon.shape[1] - 1.)
    gcp_fac = 128
    gcp_fac = numpy.maximum(gcp_fac, numpy.maximum(tie_faccell, tie_facrow))
    gcp_nrow = numpy.ceil((nrow - 1.) / gcp_fac).astype('int') + 1
    gcp_ncell = numpy.ceil((ncell - 1.) / gcp_fac).astype('int') + 1
    tie_indrow = numpy.round(numpy.linspace(0, tie_lon.shape[0] - 1,
                                            num=gcp_nrow)).astype('int')
    tie_indcell = numpy.round(numpy.linspace(0, tie_lon.shape[1] - 1,
                                             num=gcp_ncell)).astype('int')
    gcp_lon = tie_lon[tie_indrow.reshape((-1, 1)),
                      tie_indcell.reshape((1, -1))]
    gcp_lat = tie_lat[tie_indrow.reshape((-1, 1)),
                      tie_indcell.reshape((1, -1))]
    gcp_row = numpy.tile(tie_row[tie_indrow].reshape((-1, 1)) + 0.5,
                         (1, gcp_ncell))
    gcp_cell = numpy.tile(tie_cell[tie_indcell].reshape((1, -1)) + 0.5,
                          (gcp_nrow, 1))
    gcp_hei = numpy.zeros(gcp_lon.shape)
    ngcps = gcp_lon.shape

    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['gcps'] = stfmt.format_gdalgcps(gcp_lon, gcp_lat, gcp_hei,
                                                gcp_cell, gcp_row)

    syntool_stats['lon_min'] = float(numpy.min(tie_lon))
    syntool_stats['lon_max'] = float(numpy.max(tie_lon))
    syntool_stats['lat_min'] = float(numpy.min(tie_lat))
    syntool_stats['lat_max'] = float(numpy.max(tie_lat))
    return (syntool_stats, metadata, geolocation, tie_lon, tie_lat, slice_lat0,
            slice_lat1, nrow_all, ncell_all, ngcps, month)


def read_mask(infile, sltype, view, slice_lat0, slice_lat1):
    """ Read mask values from flags SLSTR file and extract mask on valid slice.
    """
    confidence_varname = 'confidence_{}{}'.format(sltype, view)
    quality_flags_filename = 'flags_{}{}.nc'.format(sltype, view)
    quality_flags_path = os.path.join(infile, quality_flags_filename)
    quality_handler = netCDF4.Dataset(quality_flags_path, 'r')
    quality_flags = quality_handler.variables[confidence_varname][slice_lat0,
                                                                  slice_lat1]
    quality_handler.close()
    return quality_flags


def read_cloud_mask(infile, sltype, view, slice_lat0, slice_lat1):
    """ Read raw cloud mask value from flags SLSTR file as the summary cloud
    flag is not reliable. """
    raw_cloud_varname = 'cloud_{}{}'.format(sltype, view)
    quality_flags_filename = 'flags_{}{}.nc'.format(sltype, view)
    quality_flags_path = os.path.join(infile, quality_flags_filename)
    quality_handler = netCDF4.Dataset(quality_flags_path, 'r')

    raw_cloud_flags = quality_handler.variables[raw_cloud_varname][slice_lat0,
                                                                   slice_lat1]
    quality_handler.close()
    return raw_cloud_flags


def read_band(infile, bandname, fieldname, slice_lat0, slice_lat1):
    """ Read band values and extract on valid slice."""
    logger.info('\tCompute contrast fort {} band'.format(bandname))
    field_filename = '{}.nc'.format(fieldname)
    data_band_path = os.path.join(infile, field_filename)
    band_handler = netCDF4.Dataset(data_band_path, 'r')
    band = band_handler.variables[fieldname][:]
    band = numpy.ma.array(band)
    band_handler.close()

    band = band[slice_lat0, slice_lat1]
    return band


def apply_default_min_max(default_minmax, bandname, vmin, vmax, syntool_stats):
    """Compute minimum and maximum values from default values if there is not
    enough valid data in the granule."""
    logger.warn('No valid values for {}'.format(bandname))
    logger.warn('Using default min/max.')
    _min = vmin
    _max = vmax
    # Use arbitrary extrema on land
    if vmin is None:
        _min = default_minmax[bandname][0]
        syntool_stats[bandname]['default_min'] = float(_min)
    if vmax is None:
        _max = default_minmax[bandname][1]
        syntool_stats[bandname]['default_max'] = float(_max)

    return _min, _max


def fromband_min_max(valid_data, bandname, vmin, vmax, syntool_stats,
                     min_percentile=0.5, max_percentile=99.99):
    """Compute histograms of valid data to retrive min and max information from
    band values."""
    _min = vmin
    _max = vmax
    if vmin is None:
        _min = numpy.percentile(valid_data, min_percentile)
        syntool_stats[bandname]['p{:.2f}'.format(min_percentile)] = float(_min)
        syntool_stats[bandname]['min'] = float(numpy.min(valid_data))
    if vmax is None:
        _max = numpy.percentile(valid_data, max_percentile)
        syntool_stats[bandname]['p{:.2f}'.format(max_percentile)] = float(_max)
        syntool_stats[bandname]['max'] = float(numpy.max(valid_data))

    return _min, _max
