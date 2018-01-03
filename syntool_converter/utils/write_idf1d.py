# encoding: utf-8

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
import math
import numpy as np
import pyproj
import logging
import argparse
from datetime import datetime
import calendar
import time
import unicodecsv
from netCDF4 import Dataset
import syntool_converter.utils.syntoolformat as stfmt

logger = logging.getLogger()
handler = logging.StreamHandler()
logger.addHandler(handler)

def write_netcdf_1d(ncfile, metadata, geolocation, band, model='trajectory',
                    ngcps=None, dgcps=None):
    """ """

    nband = len(band)
    nlin = band[0]['array'].shape[0]

    if model == 'trajectory':
        dim1name = 'time'
        unlimtime = False
    elif model == 'along_track':
        dim1name = 'time'
        unlimtime = False
    else:
        raise Exception('')
    dim1gcpname = '{}_gcp'.format(dim1name)
    if ngcps is None and dgcps is None:
       raise Exception('Need number of GCPs or pixel distance between GCPs.')
    elif ngcps is None:
        ngcps = np.ceil((nlin - 1) / float(dgcps)).astype('int') + 1
    elif dgcps is None:
        dgcps = np.ceil(nlin / (ngcps - 1.)).astype('int')
    gcplin = np.arange(ngcps).astype('int')
    gcpx = np.arange(0, nlin, 1)
    '''
    if 'geotransform' in geolocation.keys():
        srs = osr.SpatialReference()
        srs.ImportFromWkt(geolocation['projection'])
        srs4326 = osr.SpatialReference()
        srs4326.ImportFromEPSG(4326)
        if srs.IsSame(4326) == 1:
             gcplat = gcpy
             gcplon = gcpx
        else:
            proj = pyproj.Proj(srs.ExportToProj4())
            gcplon, gcplat = proj(gcpx, gcpy)
    '''
    ## Write
    dset = Dataset(ncfile, mode='w', format='NETCDF4', clobber=True)
    ## Dimensions
    if unlimtime:
        _dim1 = dset.createDimension('time', None)
    else:
        _dim1 = dset.createDimension('time', nlin)
    _dim1gcp = dset.createDimension(dim1gcpname, ngcps)
    ## Variables
    _time = dset.createVariable('time', 'f8', (dim1gcpname, ))
    _time.long_name = 'time'
    _time.standard_name = 'time'
    _time.units = 'seconds since 1970-01-01T00:00:00.000000Z'
    _time.calendar = 'standard'
    _time[:] = np.array(geolocation['time']).astype('float64')
    _latgcp = dset.createVariable('lat_gcp', 'f4', (dim1gcpname, ))
    _latgcp[:] = np.array(geolocation['geotransform'][1]).astype('float32')
    _longcp = dset.createVariable('lon_gcp', 'f4', (dim1gcpname, ))
    _longcp[:] = np.array(geolocation['geotransform'][0]).astype('float32')
    _latgcp.long_name = 'ground control points latitude'
    _latgcp.standard_name = 'latitude'
    _latgcp.units = 'degrees_north'
    _longcp.long_name = 'ground control points longitude'
    _longcp.standard_name = 'longitude'
    _longcp.units = 'degrees_east'
    _indexdim1gcp = dset.createVariable('index_{}_gcp'.format(dim1name), 'i4',
                                        (dim1gcpname, ))
    _indexdim1gcp[:] = gcplin[:].astype('int32')
    _indexdim1gcp.long_name = 'index of ground control points in {} dimension'.format(dim1name)
    _indexdim1gcp.comment = 'index goes from 0 (first pixel) to value dimension'
    for ibnd in range(nband):
        if 'name' in band[ibnd]:
            varname = band[ibnd]['name']
        elif 'description' in band[ibnd]:
            varname = band[ibnd]['description'].replace(' ', '_')
        else:
            varame = 'value_{}'.format(ibnd)
        dtype = band[ibnd]['array'].dtype
        fill_value = None
        if 'nodatavalue' in band[ibnd]:
            fill_value = dtype.type(band[ibnd]['nodatavalue'])
        _value = dset.createVariable(varname, dtype, (dim1name, ),
                                     fill_value=fill_value)
        _value[:] = band[ibnd]['array'][:]
        if 'long_name' in band[ibnd]:
            _value.long_name = band[ibnd]['long_name']
        if 'standard_name' in band[ibnd]:
            _value.standard_name = band[ibnd]['standard_name']
        if 'comment' in band[ibnd]:
            _value.comment = band[ibnd]['comment']
        if 'unittype' in band[ibnd]:
            _value.unittype = band[ibnd]['unittype']
        if 'offset' in band[ibnd]:
            _value.add_offset = np.float32(band[ibnd]['offset'])
        if 'scale' in band[ibnd]:
            _value.scale_factor = np.float32(band[ibnd]['scale'])
        if 'valid_range' in band[ibnd]:
            _value.valid_min = dtype.type(band[ibnd]['valid_range'][0])
            _value.valid_max = dtype.type(band[ibnd]['valid_range'][1])
        if 'parameter_range' in band[ibnd]:
            offset = 0
            if 'offset' in band[ibnd]:
                offset = band[ibnd]['offset']
            scale = 1
            if 'scale' in band[ibnd]:
                scale = band[ibnd]['scale']
            vran = [(p - offset) / scale for p in band[ibnd]['parameter_range']]
            vran = np.sort(np.array(vran))
            if issubclass(dtype.type, np.integer):
                vran = np.round(vran).astype(dtype)
            else:
                vran = vran.astype(dtype)
            _value.valid_min = vran[0]
            _value.valid_max = vran[1]
    ## Global attributes
    dset.idf_granule_id = metadata['name']
    dset.time_coverage_start = metadata['begin_time']
    dset.time_coverage_end = metadata['end_time']
    if model == 'trajectory':
        dset.platform_code = metadata['platform_code']
    if model == 'along_track':
        dset.cycle = metadata['cycle']
        dset.orbit = metadata['pass']
    dset.idf_subsampling_factor = np.int32(0)
    # if 'spatial_resolution' in metadata:
    #     dset.idf_spatial_resolution = metadata['spatial_resolution']
    # else:
    #    raise Exception('IDF spatial resolution is missing')
    dset.idf_spatial_resolution = np.float32(10000000)
    dset.idf_spatial_resolution_units = 'm'
    dset.close()

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument( 'input'
                       , help='Path of the CSV input file')
    parser.add_argument( 'output'
                       , help='Path of the IDF0 output file')
    parser.add_argument( '--verbose'
                       , action='store_true'
                       , default=False
                       , required=False)
    parser.add_argument( '--debug'
                       , action='store_true'
                       , default=False
                       , required=False)

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)
    if args.debug:
        logger.setLevel(logging.DEBUG)

    logger.info('test info')
    logger.debug('test debug')

    csv2idf0(args.input, args.output)
    sys.exit(0)
