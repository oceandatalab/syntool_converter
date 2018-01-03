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
import netCDF4
import logging
import syntool_converter.utils.syntoolformat as stfmt
from scipy.misc import bytescale
from scipy import interpolate
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


class OnlyNightData(Exception):
    pass


default_minmax = {'S1': (14, 150),
                  'S2': (8, 150),
                  'S3': (4, 150),
                  'S4': (None, None),
                  'S5': (None, None),
                  'S6': (None, None),
                  }


def atmospheric_correction(lut_name, infile, view, sltype, band_columns,
                           bandnames, nrow, ncell):
    '''Compute atmospheric correction from LUT '''
    # Read LUT
    sza_lut = []
    oza_lut = []
    delta_lut = []
    rho_atm = []
    oza_resol = 5
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

    # Read angles from geometry_t*.nc SLSTR
    data_angle_path = os.path.join(infile, 'geometry_t{}.nc'.format(view))
    data_angle = netCDF4.Dataset(data_angle_path, 'r')
    # Sun Zenithal Angle
    sza = data_angle.variables['solar_zenith_t{}'.format(view)][:]
    # if sza is more than 90 take values from 80 - 85 range
    sza_nan_ind = numpy.where(numpy.isnan(sza))
    sza[sza_nan_ind] = -1  # temporary, to avoid warnings with "sza > 84"
    sza[sza > 84] = 84
    sza[sza_nan_ind] = numpy.nan

    # Replace nan values by the last valid one
    ifirst = 19
    ilast = ifirst + 1
    sza[:, -ifirst:] = numpy.repeat(sza[:, -ilast],
                                    ifirst).reshape((sza.shape[0], ifirst))

    # Delta = Sun Azimuthal Angle - Observation Azimuthal Angle
    delta = (data_angle.variables['solar_azimuth_t{}'.format(view)][:]
             - data_angle.variables['sat_azimuth_t{}'.format(view)][:])
    delta = (delta + 360) % 360
    # Replace nan values by the last valid one
    delta[:, -ifirst:] = numpy.repeat(delta[:, -ilast],
                                      ifirst).reshape((delta.shape[0], ifirst))
    # Observation zenithal angle
    oza = data_angle.variables['sat_zenith_t{}'.format(view)][:]
    # Replace nan values by the last valid one
    oza[:, -ifirst:] = numpy.repeat(oza[:, -ilast],
                                    ifirst).reshape((oza.shape[0], ifirst))
    # Distance earth to sun
    data_angle.close()
    # Interpolate atmospheric correction
    nrow_geom = numpy.shape(sza)[0]
    ncell_geom = numpy.shape(sza)[1]
    # Subsampling factor across track
    ac_subsampling = int(ncell / ncell_geom)
    rho = numpy.zeros((len(band_columns), nrow, ncell))
    _points = numpy.array((sza.flatten(), oza.flatten(), delta.flatten())
                          ).transpose()
    sza_lut = numpy.unique(numpy.array(sza_lut, dtype='float32'))
    oza_lut = numpy.unique(numpy.array(oza_lut, dtype='float32'))
    delta_lut = numpy.unique(numpy.array(delta_lut, dtype='float32'))
    lut_points = (sza_lut, oza_lut, delta_lut)
    rho_atm = numpy.array(rho_atm[:][:])
    lut_shape = (sza_lut.size, oza_lut.size, delta_lut.size, 1)

    dst_range_ac = numpy.arange(ncell)
    interpolator = interpolate.RegularGridInterpolator
    for iband in range(len(band_columns)):
        lut_values = numpy.array((rho_atm[:, iband])).reshape(lut_shape)
        lut_func = interpolator(lut_points, lut_values, method='linear',
                                bounds_error=False, fill_value=None)

        _values = lut_func(_points)
        _values = _values.reshape(nrow_geom, ncell_geom)
        _values = _values * numpy.cos(sza * numpy.pi / 180.)

        src_range_ac = numpy.arange(_values.shape[1]) * ac_subsampling
        # Compensate for along track subsampling (no interpolation along track)
        if abs(nrow_geom - nrow / 2.0) > 1:
            logger.warn('Discrepency between channel and geometry dimensions')
        for irow in range(nrow_geom):
            src_values = _values[irow]
            rho[iband, irow * 2] = numpy.interp(dst_range_ac, src_range_ac,
                                                src_values)
            irownext = irow * 2 + 1
            if irownext >= nrow:
                break
            rho[iband, irownext] = numpy.interp(dst_range_ac, src_range_ac,
                                                src_values)

    # Convert reflectance into radiance TOA:
    L_toa = []
    # TODO check if there are still nan
    rho[numpy.where(numpy.isnan(rho))] = 0
    for iband in range(len(band_columns)):
        # Read Extra-terrestrial solar irradiance
        data_es_path = os.path.join(infile, '{}_quality_{}{}.nc'.format(
                                               bandnames[iband], sltype, view))
        data_es = netCDF4.Dataset(data_es_path, 'r')
        solar_flux = data_es.variables['{}_solar_irradiance_{}{}'.format(
                                          bandnames[iband], sltype, view)][:]
        # band_columns correspond to the column of the band in the LUT, which
        # is number of band + 3
        # Variation for solar flux is neglectable so we consider the mean
        es = numpy.mean(solar_flux[:])
        L_toa.append(rho[iband, :, :] * es / numpy.pi)
    data_es.close()
    return L_toa


def build_mask_rgb(channels, quality_flags, tie_lat, lat_crop):
    """"""

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

    # Mask data under sun_glint too
    off_flags = off_flags + numpy.uint32(1 << bits['sun_glint'])

    contrast_mask = (contrast_mask |
                     (numpy.bitwise_and(quality_flags, on_flags) <= 0) |
                     (numpy.bitwise_and(quality_flags, off_flags) > 0))

    # Granules containing only night/twilight data must be ignored.
    # Twilight data are only useful to improve the contrast when there is
    # at least some daily data.
    on_flags = numpy.ushort(0)
    on_flags |= 1 << bits['day']
    data_mask = (data_mask |
                 (numpy.bitwise_and(quality_flags,  on_flags) <= 0))
    if data_mask.all():
        raise OnlyNightData()

    return contrast_mask, data_mask


def get_valid_data_rgb(band, mask, max_sunglint):
    # Ignore sunglint when computing the contrast for RGB channels
    sunglint_mask = (band >= max_sunglint)
    valid_data = band.data[~mask & ~sunglint_mask]

    return valid_data


def sentinel3_slstr_rad(infile, outdir, vmin=None, vmax=None,
                        max_sunglint=150, min_percentile=2.0,
                        channels='false_rgb', write_netcdf=False,
                        lut_path=None, log_path=None, lat_crop=80.0):

    t0 = datetime.utcnow()
    # Process nadir data
    view = 'n'
    fname = 'radiance'
    sltype = 'a'
    if type(channels) is list or type(channels) is tuple:
        bandnames = channels
        product_name = 'Sentinel-3_SLSTR'
    elif 'false_rgb' == channels:
        bandnames = ('S3', 'S2', 'S1')
        product_name = 'Sentinel-3_SLSTR_false_RGB'
    else:
        raise Exception('channels must be either "false_rgb",'
                        ' or a tupple of bands')

    # convert band into column number in the LUT
    band_columns = [int(x[1:])+3-1 for x in bandnames]

    # Read coordinates and compute gcps
    (syntool_stats, metadata, geolocation, tie_lon, tie_lat, slice_lat0,
        slice_lat1, nrow_all, ncell_all, ngcps, __) = slstr.read_geometry(
                                               infile, bandnames, fname,
                                               sltype, view, product_name,
                                               vmin, vmax, log_path)

    # Compute atmospheric radiance TOA correction
    L_toa = None
    if lut_path is not None:
        t_start = datetime.utcnow()
        # Compute atmospheric correction for the whole granule
        L_toa = atmospheric_correction(lut_path, infile, view, sltype,
                                       band_columns, bandnames, nrow_all,
                                       ncell_all)
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
                                'in the geometrie_t{}.nc file'.format(view))
        t_stop = datetime.utcnow()
        syntool_stats['lut_computation'] = (t_stop - t_start).total_seconds()

    # Compute masks
    logger.info('Build masks')
    t_start = datetime.utcnow()
    quality_flags = slstr.read_mask(infile, sltype, view, slice_lat0,
                                    slice_lat1)

    try:
        contrast_mask, data_mask = build_mask_rgb(channels, quality_flags,
                                                  tie_lat, lat_crop)
    except OnlyNightData:
        logger.warn('No day data in granule.')
        sys.exit(0)
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

        # Apply atmospheric correction
        if L_toa is not None:
            band -= L_toa[band_index][:, :]

        # Mask null and negative values: they are inferior or equal to
        # atmospheric correction and should probably have been flagged as
        # clouds.
        mask_negative = (band <= 0.0)

        logger.info('\tSet contrast')
        valid_ratio_lower_threshold = 0.001  # 0.1%

        # Select valid data to compute histograms
        valid_data_mask = (band.mask | contrast_mask | mask_negative)
        valid_data = get_valid_data_rgb(band, valid_data_mask, max_sunglint)

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

        # take default min and max if min or max are too high (lake,
        # inland sea, clouds)
        _max = min(_max, max_sunglint)
        if (_min > 100):
            _min = default_minmax[bandname][0]
            _max = default_minmax[bandname][1]

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
    _min = numpy.log(_min)
    _max = numpy.log(_max)
    scale = (_max - _min) / 254.
    offset = _min
    # Construct bands
    for band_index in range(len(bandnames)):
        bandname = bandnames[band_index]
        fieldname = slstr.get_field_name(fname, bandname, sltype, view)
        band = slstr.read_band(infile, bandname, fieldname, slice_lat0,
                               slice_lat1)

        # Apply atmospheric correction
        if L_toa is not None:
            band -= L_toa[band_index][:, :]

        # Mask null and negative values: they are inferior or equal to
        # atmospheric correction and should probably have been flagged.
        mask_negative = (band <= 0.0)

        # Compute the logarithm only for radiance values that are higher than
        # the atmospheric correction.
        bnd = band.data
        bnd = numpy.log(band.data, where=(~mask_negative))

        logger.info('\tBytescaling')
        byte = bytescale(bnd, cmin=_min, cmax=_max, low=0, high=254)
        description = '{} {} (log)'.format(bandname, fname)
        if band.mask is not numpy.ma.nomask:
            byte[band.mask] = 255

        # Pixels with a radiance equal or inferior to atmospheric correction
        # are clipped to the minimal value.
        if 0 < mask_negative.size:
            byte[numpy.where(mask_negative == True)] = 0  # noqa

        # mask night data for rgb and invalid data for ir (cloud, land,
        # range value). Also mask data for extreme latitudes
        byte[data_mask] = 255

        band_range = [_vmin[band_index], _vmax[band_index]]
        bands.append({'array': byte,
                      'plot': band.data[~mask_negative],
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
        metadata['spatial_resolution'] = 500
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
