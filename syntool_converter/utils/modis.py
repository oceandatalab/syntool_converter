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

import numpy as np
import numpy.ma as ma
from pyhdf.SD import SD, SDC
from netCDF4 import Dataset
from datetime import datetime, timedelta
import ftplib
import os
import glob
import requests
import bz2
import fnmatch


# SST or OC product (http://oceancolor.gsfc.nasa.gov/cms/)
class MODISL2File(object):
    """ """
    def __init__(self, modissstfname):
        """ """
        self._fid = Dataset(modissstfname)

    def close(self):
        """ """
        self._fid.close()

    def read_sst(self):
        """ """
        fids = self._fid.groups['geophysical_data'].variables['sst']
        #att = sds.attributes()
        #scale = fidsst.att['slope']
        #offset = att['intercept']
        fillvalue = fids._FillValue #att.get('bad_value_scaled', -32767)
        sst = ma.masked_equal(fids[:, :], fillvalue).astype('float32') #* \
              #scale + offset
        return sst

    def read_chlora(self):
        """ """
        fids = self._fid.groups['geophysical_data'].variables['chlor_a']
        fillvalue = fids._FillValue #att.get('bad_value_scaled', -32767)
        chlora = ma.masked_equal(fids[:, :], fillvalue).astype('float32') 
        return chlora

    def read_lon(self):
        """ """
        return self._fid.variables('longitude')[:, :]

    def read_lat(self):
        """ """
        return self._fid.variables('latitude')[:, :]

    def read_attributes(self):
        """ """
        fmt="%Y-%m-%dT%H:%M:%S"
        attrs = {}
        #fidattrs = self._sd.attributes()
        #start_year = sdattrs['Start Year']
        #start_day = sdattrs['Start Day']
        #start_sec = sdattrs['Start Millisec'] / 1000.
        attrs['start_time'] = datetime.strptime(self._fid.time_coverage_start[:-5], fmt) #datetime(start_year, 1, 1) + \
                              #timedelta(days=start_day - 1, seconds=start_sec)
        #stop_year = sdattrs['End Year']
        #stop_day = sdattrs['End Day']
        #stop_sec = sdattrs['End Millisec'] / 1000.
        attrs['stop_time'] = datetime.strptime(self._fid.time_coverage_end[:-5], fmt)
        attrs['pass'] = self._fid.orbit_number
        #sensor = sdattrs['Sensor Name'][0:-1]
        #if 'MODISA' in sensor: # 'HMODISA' or 'MODISA'
        attrs['platform'] = self._fid.platform
        #elif 'MODIST' in sensor: # 'HMODIST' or 'MODIST'
        #    attrs['platform'] = 'Terra'
       # else:
       #     raise Exception('Unknown sensor name : {}'.format(sensor))
        return attrs


# Level-1B Calibrated Geolocated Radiances M[OY]D02
class MODIS02File(object):
    """ """
    def __init__(self, modis02fname):
        """ """
        self._sd = SD(modis02fname, SDC.READ)

    def close(self):
        """ """
        self._sd.end()

    def read_radiance(self, band):
        """ """
        sds = self._sd.select('EV_1KM_Emissive')
        att = sds.attributes()
        bands = np.array(att['band_names'].split(','), dtype='int32')
        iband = int(np.where(bands == band)[0][0])
        scale = att['radiance_scales'][iband]
        offset = att['radiance_offsets'][iband]
        fillvalue = att['_FillValue']
        rad = (ma.masked_equal(sds[iband, :, :], fillvalue).astype('float32') \
               - offset) * scale
        return rad

    def read_attributes(self):
        attrs = {}
        sdattrs = self._sd.attributes()
        coremd = sdattrs['CoreMetadata.0'].split('\n')
        for iline, line in enumerate(coremd):
            if 'OBJECT                 = RANGEENDINGDATE' in line:
                enddate = coremd[iline+2].split('"')[1]
            elif 'OBJECT                 = RANGEENDINGTIME' in line:
                endtime = coremd[iline+2].split('"')[1]
            elif 'OBJECT                 = RANGEBEGINNINGDATE' in line:
                begindate = coremd[iline+2].split('"')[1]
            elif 'OBJECT                 = RANGEBEGINNINGTIME' in line:
                begintime = coremd[iline+2].split('"')[1]
            elif 'OBJECT                 = ASSOCIATEDSENSORSHORTNAME' in line:
                attrs['sensor'] = coremd[iline+3].split('"')[1]
            elif 'OBJECT                 = ASSOCIATEDPLATFORMSHORTNAME' in line:
                attrs['platform'] = coremd[iline+3].split('"')[1]
        attrs['start_time'] = datetime.strptime(begindate+'T'+begintime,
                                                '%Y-%m-%dT%H:%M:%S.%f')
        attrs['stop_time'] = datetime.strptime(enddate+'T'+endtime,
                                               '%Y-%m-%dT%H:%M:%S.%f')
        return attrs


# Geolocation Data Set M[OY]D03
class MODIS03File(object):
    def __init__(self, modis03fname):
        """ """
        self._sd = SD(modis03fname, SDC.READ)

    def close(self):
        """ """
        self._sd.end()

    def read_lon(self):
        """ """
        return self._sd.select('Longitude')[:, :]

    def read_lat(self):
        """ """
        return self._sd.select('Latitude')[:, :]


# Cloud mask M[OY]D35_L2
class MODIS35L2File(object):
    """ """
    def __init__(self, modis35l2fname):
        """ """
        self._sd = SD(modis35l2fname, SDC.READ)

    def close(self):
        """ """
        self._sd.end()

    def read_cloudmask(self, byte=None):
        """ """
        sds = self._sd.select('Cloud_Mask')
        if byte is None:
            cloudmask = sds[:, :, :]
        else:
            cloudmask = sds[int(byte), :, :]
        return cloudmask.astype('uint8')


def modis_bright(rad, band, units):
    """ From modis_bright.pro
    """
    # Effective central wavenumber (inverse centimenters)
    cwn = [2.641775E+03, 2.505277E+03, 2.518028E+03, 2.465428E+03,
           2.235815E+03, 2.200346E+03, 0.0, 1.477967E+03,
           1.362737E+03, 1.173190E+03, 1.027715E+03, 9.080884E+02,
           8.315399E+02, 7.483394E+02, 7.308963E+02, 7.188681E+02,
           7.045367E+02]
    # Temperature correction slope (no units)
    tcs = [9.993411E-01, 9.998646E-01, 9.998584E-01, 9.998682E-01,
           9.998819E-01, 9.998845E-01, 0.0, 9.994877E-01,
           9.994918E-01, 9.995495E-01, 9.997398E-01, 9.995608E-01,
           9.997256E-01, 9.999160E-01, 9.999167E-01, 9.999191E-01,
           9.999281E-01]
    # Temperature correction intercept (Kelvin)
    tci = [4.770532E-01, 9.262664E-02, 9.757996E-02, 8.929242E-02,
           7.310901E-02, 7.060415E-02, 0.0, 2.204921E-01,
           2.046087E-01, 1.599191E-01, 8.253401E-02, 1.302699E-01,
           7.181833E-02, 1.972608E-02, 1.913568E-02, 1.817817E-02,
           1.583042E-02]
    # Compute brightness temperature
    if units == 1:
        # Radiance units are
        # Watts per square meter per steradian per micron
        result = (_bright_m(1.0e+4/cwn[band-20], rad)-tci[band-20])/tcs[band-20]
    else:
        # Radiance units are
        # milliWatts per square meter per steradian per wavenumber
        result = (_brite_m(cwn[band-20], rad)-tci[band-20])/tcs[band-20]
    return result


def _bright_m(w, r):
    """ From bright_m.pro
    """
    # Planck constant (Joule second)
    h = 6.6260755e-34
    # Speed of light in vacuum (meters per second)
    c = 2.9979246e+8
    # Boltzmann constant (Joules per Kelvin)
    k = 1.380658e-23
    # Derived constants
    c1 = 2.0*h*c*c
    c2 = (h*c)/k
    # Convert wavelength to meters
    ws = 1.0e-6*w
    # Compute brightness temperature
    return c2/(ws*np.log(c1/(1.0e+6*r*ws**5)+1.0))


def _brite_m(v, r):
    """ From brite_m.pro
    """
    # Planck constant (Joule second)
    h = 6.6260755e-34
    # Speed of light in vacuum (meters per second)
    c = 2.9979246e+8
    # Boltzmann constant (Joules per Kelvin)
    k = 1.380658e-23
    # Derived constants
    c1 = 2.0*h*c*c
    c2 = (h*c)/k
    # Convert wavenumber to inverse meters
    vs = 1.0e+2*v
    # Compute brightness temperature
    return c2*vs/np.log(c1*vs**3/(1.0e-5*r)+1.0)


def get_download_info(product_id, date):
    """ """
    year = date.strftime('%Y')
    day = date.strftime('%j')
    print(year, day, 'debug')
    if product_id.startswith('MODISAL2') or product_id.startswith('MODISTL2'):
        host = 'oceandata.sci.gsfc.nasa.gov'
        dirname = '/'.join([product_id[0:6], 'L2', year, day])
        pattern = product_id[5] + date.strftime('%Y%j%H%M%S') + \
                  '.L2_LAC_' + product_id[8:]+'.nc'
    elif product_id.startswith('MYD') or product_id.startswith('MOD'):
        host = 'ladsweb.nascom.nasa.gov'
        dirname = '/'.join(['allData', '6', product_id, year, day])
        pattern = '.'.join([product_id, 'A' + year + day, date.strftime('%H%M'),
                            '006', '*', 'hdf'])
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
    if product_id.startswith('MODISAL2') or product_id.startswith('MODISTL2'):
        url = 'https://' + host + '/cgi/getfile/' + pattern #+ '.nc'
        print 'Downloading {}'.format(url)
        urlfile = requests.get(url)
        download_path_nc = download_path #+ '.nc'
        if not os.path.exists(os.path.dirname(download_path_nc)):
            os.makedirs(os.path.dirname(download_path_nc))
        with open(download_path_nc, 'wb') as f:
            f.write(urlfile.content)
        urlfile.close()
        #print 'Extracting {}'.format(download_path_nc)
        #with open(download_path, 'wb') as f, \
        #     bz2.BZ2File(download_path_bz2, 'rb') as bz2f:
        #    f.write(bz2f.read())
            # for data in iter(lambda: bz2f.read(100 * 1024), b''):
            #     f.write(data)
        #print 'Removing {}'.format(download_path_bz2)
        #os.remove(download_path_bz2)
    elif product_id.startswith('MYD') or product_id.startswith('MOD'):
        url = os.path.join('ftp://', host, dirname, pattern)
        print 'Searching {}'.format(url)
        ftp = ftplib.FTP(host)
        ftp.login()
        ftp.cwd(dirname)
        lst = []
        ftp.retrlines('NLST', lst.append)
        flst = fnmatch.filter(lst, pattern)
        if len(flst) == 0:
            raise Exception('Could not find any {}'.format(url))
        elif len(flst) > 1:
            raise Exception('Multiple {}'.format(url))
        url = os.path.join('ftp://', host, dirname, flst[0])
        print 'Downloading {}'.format(url)
        download_path = os.path.join(download_dir, host, dirname, flst[0])
        if not os.path.exists(os.path.dirname(download_path)):
            os.makedirs(os.path.dirname(download_path))
        with open(download_path, 'wb') as f:
            ftp.retrbinary('RETR ' + flst[0], f.write)
        ftp.close()
    if not os.path.exists(download_path):
        raise Exception('{} does not exist'.format(download_path))
    return download_path
