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

from netCDF4 import Dataset
from datetime import datetime
import os
import glob
import sys
import requests
import logging

logger = logging.getLogger(__name__)

# OC product (http://oceancolor.gsfc.nasa.gov/cms/)
class VIIRSL2File(object):
    """ """
    def __init__(self, viirsl2fname):
        """ """
        self._dset = Dataset(viirsl2fname)

    def close(self):
        """ """
        self._dset.close()

    def read_chlora(self):
        """ """
        return self._dset.groups['geophysical_data'].variables['chlor_a'][:]

    def read_lon(self):
        """ """
        return self._dset.groups['navigation_data'].variables['longitude'][:]

    def read_lat(self):
        """ """
        return self._dset.groups['navigation_data'].variables['latitude'][:]

    def read_attributes(self):
        """ """
        attrs = {}
        attrs['start_time'] = datetime.strptime(self._dset.time_coverage_start,
                                                '%Y-%m-%dT%H:%M:%S.%fZ')
        attrs['stop_time'] = datetime.strptime(self._dset.time_coverage_end,
                                               '%Y-%m-%dT%H:%M:%S.%fZ')
        attrs['pass'] = self._dset.startDirection
        attrs['platform'] = self._dset.platform
        return attrs


def get_download_info(product_id, date):
    """ """
    year = date.strftime('%Y')
    day = date.strftime('%j')
    if product_id == 'VIIRSL2OC':
        host = 'oceandata.sci.gsfc.nasa.gov'
        dirname = '/'.join(['VIIRS', 'L2', year, day])
        pattern = 'V' + date.strftime('%Y%j%H%M%S') + '.L2_SNPP_OC.nc'
    else:
        raise Exception('Unknown product_id')
    return (host, dirname, pattern)


def search_and_download(product_id, date, download_dir):
    """ """
    host, dirname, pattern = get_download_info(product_id, date)
    # Search
    download_path = os.path.join(download_dir, host, dirname, pattern)
    print 'Searching {}'.format(download_path)
    found = glob.glob(download_path)
    if len(found) > 0:
        found.sort()
        return found[-1]
    # Download
    if product_id == 'VIIRSL2OC':
        url = 'https://' + host + '/cgi/getfile/' + pattern
        print 'Downloading {}'.format(url)
        r = None
        try:
            r = requests.get(url)
        except requests.ConnectionError:
            # Handle "unreachable server" errors
            _, e, _ = sys.exc_info()
            logger.error('Remote server "{}" is unreachable'.format(host))
            logger.debug(e)
            raise
        if not r.ok:
            # Handle HTTP-related errors
            logger.error(r.text)
            raise
        if not os.path.exists(os.path.dirname(download_path)):
            os.makedirs(os.path.dirname(download_path))
        with open(download_path, 'wb') as f:
            f.write(r.content)
    if not os.path.exists(download_path):
        raise Exception('{} does not exist'.format(download_path))
    return download_path
