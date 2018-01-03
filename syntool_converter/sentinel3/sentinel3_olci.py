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
import scipy
import netCDF4
import logging
import syntool_converter.utils.syntoolformat as stfmt
from scipy.misc import bytescale
from datetime import datetime

logger = logging.getLogger(__name__)

# S-2 : 665, 560, 490
# red : 650 nm
# green : 510 nm
# blue : 475 nm
# https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci/overview/heritage
# https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-3-olci/level-1/fr-or-rr-toa-radiances
# Oa1   400 	15
# Oa2 	412.5 	10
# Oa3 	442.5 	10
# Oa4 	490 	10   BLUE
# Oa5 	510 	10
# Oa6 	560 	10   GREEN
# Oa7 	620 	10
# Oa8 	665 	10   RED
# Oa9 	673.75 	7.5
# Oa10 	681.25 	7.5
# Oa11 	708.75 	10
# Oa12 	753.75 	7.5
# Oa13 	761.25 	2.5
# Oa14 	764.375 3.75
# Oa15 	767.5 	2.5
# Oa16 	778.75 	15
# Oa17 	865 	20 SUN GLITTER
# Oa18 	885 	10
# Oa19 	900 	10
# Oa20 	940 	20
# Oa21 	1 020 	40

default_minmax = {'Oa01': (None, None),
                  'Oa02': (None, None),
                  'Oa03': (None, None),
                  'Oa04': (13.966, 95.000),
                  'Oa05': (None, None),
                  'Oa06': (8.213, 95.000),
                  'Oa07': (None, None),
                  'Oa08': (None, None),
                  'Oa09': (3.346, 95.000),
                  'Oa10': (None, None),
                  'Oa11': (None, None),
                  'Oa12': (None, None),
                  'Oa13': (None, None),
                  'Oa14': (None, None),
                  'Oa15': (None, None),
                  'Oa16': (None, None),
                  'Oa17': (3., 150.0),
                  'Oa18': (None, None),
                  'Oa19': (None, None),
                  'Oa20': (None, None),
                  'Oa21': (None, None)}


def atmospheric_correction(lut_name, infile, band_columns, nrow, ncell):
    '''Compute atmospheric correction from LUT '''
    # Read LUT
    sza_lut = []
    oza_lut = []
    delta_lut = []
    rho_atm = []
    oza_resol = 5
    # sza_resol = 5
    delta_resol = 10
    ntot_delta = 360 / delta_resol + 1
    ntot_oza = 90 / oza_resol + 1
    with open(lut_name, 'r') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            split_line = line.split(' ')
            sza_lut.append(split_line[0])
            oza_lut.append(split_line[1])
            delta_lut.append(split_line[2])
            rho_line = []
            for icolband in band_columns:
                rho_line.append(float(split_line[icolband]))
            rho_atm.append(rho_line)

    # Check the size of the LUT
    lut_size = len(rho_atm)
    lut_size_factor = ntot_oza * ntot_delta
    if 0 != (lut_size % lut_size_factor):
        raise Exception('The size of the LUT does not match the hardcoded '
                        'parameters: the number of lines per sza value should '
                        'be {}'.format(lut_size_factor))

    # Read angles from tie_geometrie from OLCI
    data_angle_path = os.path.join(infile, 'tie_geometries.nc')
    data_angle = netCDF4.Dataset(data_angle_path, 'r')
    # Sun Zenithal Angle
    sza = data_angle.variables['SZA'][:]
    # Delta = Sun Azimuthal Angle - Observation Azimuthal Angle
    delta = data_angle.variables['SAA'][:] - data_angle.variables['OAA'][:]
    delta = (delta + 360) % 360
    # Observation zenithal angle
    oza = data_angle.variables['OZA'][:]
    # Subsampling factor across track
    ac_subsampling = data_angle.ac_subsampling_factor
    data_angle.close()
    # Interpolate atmospheric correction
    nrow_geom = numpy.shape(sza)[0]
    ncell_geom = numpy.shape(sza)[1]
    if nrow_geom != nrow:
        raise Exception('Wrong along track dimension in tie_geometries.nc')
    rho = numpy.zeros((len(band_columns), nrow, ncell))
    _points = numpy.array((sza.flatten(), oza.flatten(), delta.flatten())
                          ).transpose()
    sza_lut = numpy.unique(numpy.array(sza_lut, dtype='float32'))
    oza_lut = numpy.unique(numpy.array(oza_lut, dtype='float32'))
    delta_lut = numpy.unique(numpy.array(delta_lut, dtype='float32'))
    lut_points = (sza_lut, oza_lut, delta_lut)
    rho_atm = numpy.array(rho_atm[:][:])
    lut_shape = (sza_lut.size, oza_lut.size, delta_lut.size, 1)

    dst_range = numpy.arange(ncell)
    interpolator = scipy.interpolate.RegularGridInterpolator
    if nrow != nrow_geom:
        logger.warn('Discrepency between channel and geometry dimensions')
    for iband in range(len(band_columns)):
        lut_values = numpy.array((rho_atm[:, iband])).reshape(lut_shape)
        lut_func = interpolator(lut_points, lut_values, method='linear',
                                bounds_error=False, fill_value=None)

        _values = lut_func(_points)
        _values = _values.reshape(nrow_geom, ncell_geom)
        _values = _values * numpy.cos(sza * numpy.pi / 180.)

        src_range = numpy.arange(_values.shape[1]) * ac_subsampling
        for irow in range(nrow_geom):
            src_values = _values[irow]
            if irow >= nrow:
                break
            rho[iband, irow] = numpy.interp(dst_range, src_range, src_values)

    # Convert reflectance into radiance TOA:
    # Read Extra-terrestrial solar irradiance
    data_es_path = os.path.join(infile, 'instrument_data.nc')
    data_es = netCDF4.Dataset(data_es_path, 'r')
    solar_flux = data_es.variables['solar_flux'][:]
    data_es.close()
    # 10**(-3) variation for solar flux is neglectable so we consider the mean
    L_toa = []
    for iband in range(len(band_columns)):
        # band_columns correspond to the column of the band in the LUT, which
        # is number of band + 3
        es = numpy.mean(solar_flux[band_columns[iband]-3, :])
        L_toa.append(rho[iband, :, :] * es / numpy.pi)
    return L_toa


def sentinel3_olci(infile, outdir, vmin=None, vmax=None,
                   slope_threshold=-0.00001, channels='nir',
                   write_netcdf=False, lut_path=None, log_path=None,
                   lat_crop=85.0):
    """"""
    t0 = datetime.utcnow()

    if type(channels) is list or type(channels) is tuple:
        bandnames = channels
        product_name = 'Sentinel-3_OLCI'
    elif 'nir' == channels:
        bandnames = ('Oa17',)
        product_name = 'Sentinel-3_OLCI_NIR'
    elif 'true_rgb' == channels:
        bandnames = ('Oa09', 'Oa06', 'Oa04')
        product_name = 'Sentinel-3_OLCI_true_RGB'
    elif 'false_rgb' == channels:
        bandnames = ('Oa17', 'Oa06', 'Oa04')
        product_name = 'Sentinel-3_OLCI_false_RGB'
    else:
        raise Exception('channels must be either "nir", "true_rgb" '
                        'or "false_rgb"')

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

    # convert band into column number in the LUT
    band_columns = [int(x[2:])+3-1 for x in bandnames]

    if vmin is None:
        vmin = [None] * len(bandnames)
    if vmax is None:
        vmax = [None] * len(bandnames)

    full_path = os.path.normpath(infile)
    file_path = os.path.basename(full_path)
    file_name, _ = os.path.splitext(file_path)

    geo_path = os.path.join(infile, 'geo_coordinates.nc')
    time_path = os.path.join(infile, 'time_coordinates.nc')
    quality_path = os.path.join(infile, 'qualityFlags.nc')

    # Extract geo coordinates information
    geo_handler = netCDF4.Dataset(geo_path, 'r')
    nrow = geo_handler.dimensions['rows'].size
    nrow_all = nrow
    ncell = geo_handler.dimensions['columns'].size
    ncell_all = ncell
    lon = geo_handler.variables['longitude'][:]
    tie_lon = numpy.ma.array(lon)
    lat = geo_handler.variables['latitude'][:]
    tie_lat = numpy.ma.array(lat)
    geo_handler.close()

    # Handle longitude continuity
    dlon = lon[1:, :] - lon[:-1, :]
    if 180.0 <= numpy.max(numpy.abs(dlon)):
        lon[lon < 0.0] = lon[lon < 0.0] + 360.0

    # Extract time coordinates information
    time_handler = netCDF4.Dataset(time_path, 'r')
    start_timestamp = time_handler.variables['time_stamp'][0]
    end_timestamp = time_handler.variables['time_stamp'][-1]
    timestamp_units = time_handler.variables['time_stamp'].units
    time_handler.close()

    # Format time information
    start_time = netCDF4.num2date(start_timestamp, timestamp_units)
    end_time = netCDF4.num2date(end_timestamp, timestamp_units)
    (dtime, time_range) = stfmt.format_time_and_range(start_time, end_time,
                                                      units='ms')

    parameters = ['{} TOA radiance'.format(bnd) for bnd in bandnames]
    metadata = {}
    metadata['product_name'] = product_name
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
    metadata['sensor_type'] = 'medium-resolution imaging spectrometer'
    metadata['sensor_name'] = 'OLCI'
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

    geolocation = {}
    geolocation['projection'] = stfmt.format_gdalprojection()
    geolocation['gcps'] = stfmt.format_gdalgcps(gcp_lon, gcp_lat, gcp_hei,
                                                gcp_cell, gcp_row)

    syntool_stats['lon_min'] = float(numpy.min(tie_lon))
    syntool_stats['lon_max'] = float(numpy.max(tie_lon))
    syntool_stats['lat_min'] = float(numpy.min(tie_lat))
    syntool_stats['lat_max'] = float(numpy.max(tie_lat))

    # Compute atmospheric radiance TOA correction
    L_toa = None
    if lut_path is not None:
        t_start = datetime.utcnow()
        # Compute atmospheric correction for the whole granule
        L_toa = atmospheric_correction(lut_path, infile, band_columns,
                                       nrow_all, ncell_all)
        # Extract the slices of atmospheric correction that match the ones of
        # the bands after the removal of the high latitude columns.
        for iband in range(len(L_toa)):
            L_toa[iband] = L_toa[iband][slice_lat0, slice_lat1]
            # An erroneous LUT could produce nan values for the atmospheric
            # correction. This is not acceptable.
            if numpy.any(~numpy.isfinite(L_toa[iband])):
                raise Exception('Infinite or nan value found in the '
                                'atmospheric correction, please check that the'
                                'LUT has valid values for the angles contained'
                                'in the instrument_data.nc file')
        t_stop = datetime.utcnow()
        syntool_stats['lut_computation'] = (t_stop - t_start).total_seconds()

    logger.info('Construct bands')
    t_start = datetime.utcnow()
    bands = []
    dset = netCDF4.Dataset(quality_path, 'r')
    quality_flags = dset.variables['quality_flags'][slice_lat0, slice_lat1]
    dset.close()

    bits = {'saturated@Oa21': 0,
            'saturated@Oa20': 1,
            'saturated@Oa19': 2,
            'saturated@Oa18': 3,
            'saturated@Oa17': 4,
            'saturated@Oa16': 5,
            'saturated@Oa15': 6,
            'saturated@Oa14': 7,
            'saturated@Oa13': 8,
            'saturated@Oa12': 9,
            'saturated@Oa11': 10,
            'saturated@Oa10': 11,
            'saturated@Oa09': 12,
            'saturated@Oa08': 13,
            'saturated@Oa07': 14,
            'saturated@Oa06': 15,
            'saturated@Oa05': 16,
            'saturated@Oa04': 17,
            'saturated@Oa03': 18,
            'saturated@Oa02': 19,
            'saturated@Oa01': 20,
            'dubious': 21,
            'sun-glint_risk': 22,
            'duplicated': 23,
            'cosmetic': 24,
            'invalid': 25,
            'straylight_risk': 26,
            'bright': 27,
            'tidal_region': 28,
            'fresh_inland_water': 29,
            'coastline': 30,
            'land': 31}

    off_flags = numpy.uint32(0)
    off_flags = off_flags + numpy.uint32(1 << bits['dubious'])
    off_flags = off_flags + numpy.uint32(1 << bits['invalid'])
    off_flags = off_flags + numpy.uint32(1 << bits['straylight_risk'])
    off_flags = off_flags + numpy.uint32(1 << bits['bright'])
    off_flags = off_flags + numpy.uint32(1 << bits['tidal_region'])
    # off_flags = off_flags + numpy.uint32(1 << bits['fresh_inland_water'])
    off_flags = off_flags + numpy.uint32(1 << bits['coastline'])
    off_flags = off_flags + numpy.uint32(1 << bits['land'])

    # Filter out values where any of the bands is flagged as saturated
    for bandname in bandnames:
        bit_name = 'saturated@{}'.format(bandname)
        off_flags = off_flags + numpy.uint32(1 << bits[bit_name])

    lat_mask = (numpy.abs(tie_lat) > lat_crop)
    data_mask = numpy.zeros(quality_flags.shape, dtype='bool')
    data_mask = (data_mask | lat_mask)

    if numpy.all(data_mask):
        logger.warn('No data to extract.')
        sys.exit(0)

    contrast_mask = numpy.zeros(quality_flags.shape, dtype='bool')
    contrast_mask = (contrast_mask | lat_mask)
    contrast_mask = (contrast_mask |
                     (numpy.bitwise_and(quality_flags, off_flags) > 0))
    t_stop = datetime.utcnow()
    syntool_stats['mask_computation'] = (t_stop - t_start).total_seconds()

    t_start = datetime.utcnow()
    _vmin = list(vmin)
    _vmax = list(vmax)

    for band_index in range(len(bandnames)):
        bandname = bandnames[band_index]
        logger.info('\tConstruct {} band'.format(bandname))
        fieldname = '{}_radiance'.format(bandname)
        file_path = os.path.join(infile, '{}.nc'.format(fieldname))
        f_handler = netCDF4.Dataset(file_path, 'r')
        band = f_handler.variables[fieldname][:]
        band = numpy.ma.array(band)
        f_handler.close()

        # Apply atmospheric correction
        band = band[slice_lat0, slice_lat1]
        if L_toa is not None:
            band -= L_toa[band_index][:, :]

        # Mask null and negative values: they are inferior or equal to
        # atmospheric correction and should probably have be flagged as clouds.
        mask_negative = (band <= 0.0)

        logger.info('\tSet contrast')
        valid_ratio_lower_threshold = 0.001  # 0.1%
        valid_data_mask = (band.mask | contrast_mask | mask_negative)
        valid_data = band.data[~valid_data_mask]
        valid_ratio = float(valid_data.size) / float(band.data.size)
        syntool_stats[bandname]['valid_ratio'] = valid_ratio
        if valid_ratio_lower_threshold >= valid_ratio:
            logger.warn('No valid values for {}'.format(bandname))
            logger.warn('Using default min/max.')

            # Use arbitrary extrema on land
            if _vmin[band_index] is None:
                _min = default_minmax[bandname][0]
                _vmin[band_index] = _min
                syntool_stats[bandname]['default_min'] = float(_min)
            if _vmax[band_index] is None:
                _max = default_minmax[bandname][1]
                _vmax[band_index] = _max
                syntool_stats[bandname]['default_max'] = float(_max)
        else:
            # TODO: add clipping for NIR
            # _min = numpy.clip(_min, 1.5, 10.0)
            # _max = numpy.clip(_max, 30.0, 60.0)
            if _vmin[band_index] is None:
                _min = numpy.percentile(valid_data, .5)
                _vmin[band_index] = _min
                syntool_stats[bandname]['p0050'] = float(_min)
                syntool_stats[bandname]['min'] = float(numpy.min(valid_data))
            if _vmax[band_index] is None:
                _max = numpy.percentile(valid_data, 99.99)
                _vmax[band_index] = _max
                p99 = numpy.percentile(valid_data, 99.0)
                syntool_stats[bandname]['p9900'] = float(p99)
                syntool_stats[bandname]['p9999'] = float(_max)
                syntool_stats[bandname]['max'] = float(numpy.max(valid_data))
        logger.info('\tContrast : vmin={} / vmax={}'.format(_vmin[band_index],
                                                            _vmax[band_index]))

    min_values = [_vmin[band_index] for band_index in range(len(bandnames))]
    max_values = [_vmax[band_index] for band_index in range(len(bandnames))]

    t_stop = datetime.utcnow()
    syntool_stats['minmax_computation'] = (t_stop - t_start).total_seconds()
    syntool_stats['final_min'] = float(numpy.min(min_values))
    syntool_stats['final_max'] = float(numpy.max(max_values))

    _min = numpy.log(numpy.min(min_values))
    _max = numpy.log(numpy.max(max_values))
    scale = (_max - _min) / 254.
    offset = _min
    for band_index in range(len(bandnames)):
        bandname = bandnames[band_index]
        logger.info('\tConstruct {} band'.format(bandname))
        fieldname = '{}_radiance'.format(bandname)
        file_path = os.path.join(infile, '{}.nc'.format(fieldname))
        f_handler = netCDF4.Dataset(file_path, 'r')
        band = f_handler.variables[fieldname][:]
        band = numpy.ma.array(band)
        f_handler.close()

        # Apply atmospheric correction
        band = band[slice_lat0, slice_lat1]
        if L_toa is not None:
            band -= L_toa[band_index][:, :]

        # Mask null and negative values: they are inferior or equal to
        # atmospheric correction and should probably have be flagged as clouds.
        mask_negative = (band <= 0.0)

        # Compute the logarithm only for radiance values that are higher than
        # the atmospheric correction.
        bnd = numpy.log(band.data, where=(~mask_negative))

        logger.info('\tBytescaling')
        byte = bytescale(bnd, cmin=_min, cmax=_max, low=0, high=254)
        description = '{} radiance (log)'.format(bandname)
        if band.mask is not numpy.ma.nomask:
            byte[band.mask] = 255

        # mask data for extreme latitudes
        byte[data_mask] = 255

        # Pixels with a radiance equal or inferior to atmospheric correction
        # are clipped to the minimal value.
        if 0 < mask_negative.size:
            byte[numpy.where(mask_negative == True)] = 0  # noqa

        band_range = [_vmin[band_index], _vmax[band_index]]
        bands.append({'array': byte,
                      'scale': scale,
                      'offset': offset,
                      'description': description,
                      'unittype': '',
                      'nodatavalue': 255,
                      'parameter_range': band_range})
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
        metadata['spatial_resolution'] = 300
        ncfile = stfmt.format_ncfilename(outdir, metadata, create_dir=True)
        stfmt.write_netcdf(ncfile, metadata, geolocation, bands, 'swath',
                           ngcps=gcp_lon.shape)
    else:
        logger.info('Write geotiff')
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, bands)

    logger.info(datetime.utcnow() - t0)
    syntool_stats['total_time'] = (datetime.utcnow() - t0).total_seconds()
    if log_path is not None:
        import json
        stats_path = os.path.join(log_path, '{}.json'.format(file_name))
        with open(stats_path, 'w') as f:
            json.dump(syntool_stats, f)
