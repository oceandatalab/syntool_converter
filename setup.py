# -*- encoding: utf-8 -*-

"""
syntool_converter: conversion of raw data files to Syntool GeoTiff format.

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
from setuptools import setup
import subprocess

package_dir = os.path.dirname(__file__)
version_path = os.path.join(package_dir, 'VERSION.txt')

major_version = '0.1'
if os.path.exists('.git') and os.path.isdir('.git'):
    commits = subprocess.check_output([ '/usr/bin/git'
                                      , 'rev-list'
                                      , 'HEAD'
                                      , '--count']).decode('utf-8').strip()
    with open(version_path, 'w') as f:
        f.write('{}.{}\n'.format(major_version, commits))

with open(version_path, 'r') as f:
    version = f.read()

setup(
    zip_safe=False,
    name='syntool_converter',
    version=version,
    author=', '.join(('Gilles Guitton <gilles.guitton@oceandatalab.com>',
                     'Lucile Gaultier <lucile.gaultier@oceandatalab.com>',
                     'Sylvain Herlédan <sylvain.herledan@oceandatalab.com>')),
    author_email='syntool@oceandatalab.com',
    packages=[ 'syntool_converter'
             , 'syntool_converter.bathy'
             , 'syntool_converter.current'
             , 'syntool_converter.duacs'
             , 'syntool_converter.model_wave'
             , 'syntool_converter.model_wind'
             , 'syntool_converter.modis'
             , 'syntool_converter.sar'
             , 'syntool_converter.sea_ice'
             , 'syntool_converter.scat'
             , 'syntool_converter.sentinel2'
             , 'syntool_converter.sentinel3'
             , 'syntool_converter.smos'
             , 'syntool_converter.smosstorm'
             , 'syntool_converter.sst'
             , 'syntool_converter.utils'
             , 'syntool_converter.viirs'
             ],
    scripts=[
        'bin/syntool-converter',
        'bin/syntool-colorbar',
    ],
    url='https://git.oceandatalab.com/syntool_odl/syntool_converter',
    license='AGPLv3',
    description='Conversion of raw data files to Syntool GeoTiff format.',
    long_description=open('README.txt').read(),
    install_requires=[ 'gdal'
                     , 'Pillow'
                     , 'numpy'
                     , 'scipy'
                     , 'matplotlib'
                     , 'pyproj'
                     , 'pyresample'
                     , 'netCDF4'
                     , 'python-hdf4'
                     , 'ceraux'
                     , 'cerbere'
                     , 'sar'
                     , 'requests[security]'
    ],
    package_data={ 'syntool_converter.utils': ['palette/*.*']
                 , 'syntool_converter.model_wind': ['ecmwf_0125_land_mask.nc']
                 , 'syntool_converter': ['share/colorbars/config.json',
                                         'share/range_files/*.txt']
    },
)
