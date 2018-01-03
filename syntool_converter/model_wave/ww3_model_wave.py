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

from cerbere.mapper.ncfile import NCFile
import numpy as np
from netCDF4 import num2date, date2num
from datetime import datetime
import os
import syntool_converter.utils.syntoolformat as stfmt
import osr


# def show(infile, name='phs'):
#     """
#     """
#     import matplotlib.pyplot as plt
#     if name == 'phs':
#         varnames = ['phs0', 'phs1', 'phs2', 'phs3']
#     elif name == 'ptp':
#         varnames = ['ptp0', 'ptp1', 'ptp2', 'ptp3']
#     elif name == 'pdir':
#         varnames = ['pdir0', 'pdir1', 'pdir2', 'pdir3']
#     elif name == 'pws':
#         varnames = ['pws0', 'pws1', 'pws2', 'pws3']
#     fig = plt.figure()
#     ww3 = NCFile(infile)
#     vars = [ww3.read_values(varname)[0, :, :] for varname in varnames]
#     vmax = np.array(vars).max()
#     for i, var in enumerate(vars):
#         plt.subplot(2, 2, i+1)
#         plt.imshow(var, origin='lower', vmax=vmax)
#         plt.colorbar()
#         plt.title(varnames[i])
#     plt.show()


# def arctic_grid():
#     """
#     """
#     #fname = '/local/home/data/ancillary/ww3/HINDCAST/ARCTIC/2014_ECMWF_CERSAT/hs/ww3.201409_hs.nc'
#     fname = '/home/cercache/project/ww3/public/HINDCAST/ARCTIC/2014_ECMWF_new/partitions/ww3.arc06m.201404_phs0.nc'
#     ncfile = NCFile(fname)
#     lon = ncfile.read_values('lon')
#     lat = ncfile.read_values('lat')
#     srs = osr.SpatialReference()
#     srs.ImportFromEPSG(3411)
#     import pyproj
#     proj = pyproj.Proj(srs.ExportToProj4())
#     x, y = proj(lon, lat)
#     x0, x1 = x[:, 0].mean(), x[:, -1].mean()
#     dx = (x1-x0)/(lon.shape[1]-1)
#     y0, y1 = y[0, :].mean(), y[-1, :].mean()
#     dy = (y1-y0)/(lon.shape[0]-1)
#     print x0, x1, dx, x0-dx/2, x1+dx/2
#     print y0, y1, dy, y0-dy/2, y1+dy/2
#     import pdb ; pdb.set_trace()


def split_ww3_fname(infile, varname):
    """
    """
    import re
    dirname = os.path.dirname(os.path.dirname(infile))
    basename = os.path.basename(infile).rsplit('_', 1)[0] + '_{}.nc'
    if re.match(r'((phs)|(ptp)|(pdir))[0-9]', varname) != None:
        fname = os.path.join(dirname, 'partitions', basename.format(varname))
    else:
        fname = os.path.join(dirname, varname, basename.format(varname))
    return fname


def other_split_ww3_fname(infile, varname):
    """
    """
    dirname = os.path.dirname(infile)
    basename = os.path.basename(infile).split('_')[0] + '_{}.nc'
    fname = os.path.join(dirname, basename.format(varname))
    return fname


def ww3_model_wave(infile, outdir, date=None, max_forecast_hours=None,
                   vmin=0., vmax=25.4, vmin_pal=0., vmax_pal=10., v2=False):
    """
    """
    # Read/Process data
    print 'Read/Process data'
    ncfile = NCFile(infile)
    ww3 = {}
    ww3['time'] = ncfile.read_values('time')
    ww3['time_units'] = ncfile.read_field('time').units
    ww3uniqtime = ww3['time'].size == 1
    if ww3uniqtime: # assume 3h
        ww3['deltatime'] = 180
    else:
        t01 = num2date(np.array(ww3['time'][0:2]), ww3['time_units'])
        ww3['deltatime'] = np.round((t01[1] - t01[0]).total_seconds() / 60.).astype('int')
    if date != None:
        tind = np.where(ww3['time'] == date2num(date, ww3['time_units']))[0]
        if tind.size == 1:
            tsl = slice(tind[0], tind[0]+1)
            ww3['time'] = ww3['time'][tsl]
        else:
            raise Exception('Date not found in WW3 file.')
    else:
        tsl = slice(0, ww3['time'].size)
    if max_forecast_hours is not None:
        if ww3['time'].size != 1:
            raise Exception('max_forecast_hours option works with only 1 time.')
        ww3time = num2date(ww3['time'][0], ww3['time_units'])
        if 'date_cycle' not in ncfile.read_global_attributes():
            raise Exception('max_forecast_hours option works with date_cycle attribute.')
        cycletime_str = ncfile.read_global_attribute('date_cycle')
        if 10 == len(cycletime_str):
            cycletime_format = '%Y%m%d%H'
        elif 20 == len(cycletime_str) and cycletime_str.endswith('Z'):
            cycletime_format = '%Y-%m-%dT%H:%M:%SZ'
        else:
            raise Exception('Cycletime format is not supported: {}'.format(cycletime_str))

        cycletime = datetime.strptime(cycletime_str, cycletime_format)
        forecast_hours = (ww3time - cycletime).total_seconds() / 3600.
        if forecast_hours > max_forecast_hours:
            raise Exception('Exceeds max_forecast_hours.')
    slices = [tsl,
              slice(0, ncfile.get_dimsize('latitude')),
              slice(0, ncfile.get_dimsize('longitude'))]
    # slices = [tsl,
    #           slice(ncfile.get_dimsize('latitude')-1, None, -1),
    #           slice(0, 2*ncfile.get_dimsize('longitude'))]
    ncfieldnames = ncfile.get_fieldnames()
    fieldnames = ['hs', 'phs0', 'phs1', 'phs2', 'phs3',
                  'ptp0', 'ptp1', 'ptp2', 'ptp3',
                  'pdir0', 'pdir1', 'pdir2', 'pdir3']
    ww3['source'] = [infile]
    for fieldname in fieldnames:
        if fieldname in ncfieldnames:
            ww3[fieldname] = ncfile.read_values(fieldname, slices=slices)[:, ::-1, :]
        else:
            infile2 = split_ww3_fname(infile, fieldname)
            if not os.path.exists(infile2):
                infile2 = other_split_ww3_fname(infile, fieldname)
            if os.path.exists(infile2):
                ncfile2 = NCFile(infile2)
                ww3[fieldname] = ncfile2.read_values(fieldname, slices=slices)[:, ::-1, :]
                ncfile2.close()
                ww3['source'].append(infile2)
    ww3['npart'] = 0
    for indp in range(4):
        pin = [name+str(indp) in ww3 for name in ['phs', 'ptp', 'pdir']]
        if all(pin):
            ww3['npart'] += 1
        else:
            break
    if ww3['npart'] == 0:
        raise Exception('Could not find all partition variables.')
    ww3['area'] = ncfile.read_global_attribute('area')
    if 'global' in ww3['area'].lower():
        ww3['lon'] = ncfile.read_values('lon')
        ww3['lat'] = ncfile.read_values('lat')[::-1]
        ww3['lon_res'] = float(ncfile.read_global_attribute('longitude_resolution'))
        ww3['lat_res'] = float(ncfile.read_global_attribute('latitude_resolution'))
    elif ww3['area'] == 'ARCTIC-12km':
        ww3['lon'] = ncfile.read_values('lon')[::-1, :]
        ww3['lat'] = ncfile.read_values('lat')[::-1, :]
    else:
        raise Exception('Not implemented : area = "{}"'.format(ww3['area']))
    if 'run_time' in ncfile.read_global_attributes():
        run_time = ncfile.read_global_attribute('run_time')
        ww3['rundtime'] = datetime.strptime(run_time, '%Y-%m-%dT%H:%M:%SZ')
    ncfile.close()

    # Construct metadata/geolocation/band(s)
    print 'Construct metadata/geolocation/band(s)'
    metadata = {}
    #metadata['product_name'] =
    #metadata['name'] =
    #metadata['datetime'] =
    metadata['time_range'] = ['-{:d}m'.format(ww3['deltatime'] / 2),
                              '+{:d}m'.format(ww3['deltatime'] / 2)]
    metadata['source_URI'] = ww3['source']
    metadata['source_provider'] = ['SHOM', 'Ifremer']
    metadata['processing_center'] = ''
    metadata['conversion_software'] = 'Syntool'
    metadata['conversion_version'] = '0.0.0'
    metadata['conversion_datetime'] = stfmt.format_time(datetime.utcnow())
    #metadata['parameter'] =
    metadata['type'] = 'model'
    #metadata['model_longitude_resolution'] =
    #metadata['model_latitude_resolution'] =
    if 'rundtime' in ww3:
        metadata['model_analysis_datetime'] = stfmt.format_time(ww3['rundtime'])
    geolocation = {}
    if 'global' in ww3['area'].lower():
        metadata['model_longitude_resolution'] = ww3['lon_res']
        metadata['model_latitude_resolution'] = ww3['lat_res']
        geolocation['projection'] = stfmt.format_gdalprojection()
        lon0, dlon = ww3['lon'][0], ww3['lon'][1] - ww3['lon'][0]
        lat0, dlat = ww3['lat'][0], ww3['lat'][1] - ww3['lat'][0]
        geolocation['geotransform'] = [lon0 - dlon / 2., dlon, 0,
                                       lat0 - dlat / 2., 0, dlat]
        ww3['grid'] = 'GLOBAL'
    elif ww3['area'] == 'ARCTIC-12km':
        import pyproj
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3411)
        proj = pyproj.Proj(srs.ExportToProj4())
        x01, dummy = proj(ww3['lon'][:, [0, -1]], ww3['lat'][:, [0, -1]])
        x0, x1 = x01.mean(axis=0)
        dx = (x1 - x0) / (ww3['lon'].shape[1] - 1)
        dummy, y01 = proj(ww3['lon'][[0, -1], :], ww3['lat'][[0, -1], :])
        y0, y1 = y01.mean(axis=1)
        dy = (y1 - y0) / (ww3['lon'].shape[0] - 1)
        geolocation['projection'] = srs.ExportToWkt()
        geolocation['geotransform'] = [x0 - dx / 2., dx, 0, y0 - dy / 2., 0, dy]
        # geolocation['geotransform'] = [-2600051.73564, 12500.2285676, 0,
        #                                2787547.79214, 0, -12500.2262608]
        ww3['grid'] = 'ARCTIC'
    else:
        raise Exception('Not implemented : area = "{}"'.format(ww3['area']))

    # Loop on time
    for itime in range(ww3['time'].size):
        dtime = num2date(ww3['time'][itime], ww3['time_units'])
        metadata['datetime'] = stfmt.format_time(dtime)
        if len(ww3['source']) == 1:
            basename = os.path.splitext(os.path.basename(infile))[0]
            if not ww3uniqtime:
                _date = dtime.strftime('%Y%m%d')
                _datetime = dtime.strftime('%Y%m%dT%H')
                if _date in basename and _datetime not in basename:
                    basename = basename.replace(_date, dtime.strftime('%Y%m%dT%HZ'))
                else:
                    raise Exception
        else:
            basename = 'WW3-' + ww3['grid'] + '-' + dtime.strftime('%Y%m%dT%HZ')
            if 'hindcast' in infile.lower():
                basename = 'HINDCAST_' + basename

        ### Total HS ###
        # Update metadata
        metadata['product_name'] = 'WW3_model_wave_hs'
        if v2 == True:
            metadata['product_name'] += '_v2'
        metadata['name'] = basename + '-hs'
        metadata['parameter'] = 'wave significant height'
        # Make band
        band = []
        hs = ww3['hs'][itime, :, :]
        offset, scale = vmin, (vmax - vmin) / 254.0
        np.clip(hs, vmin, vmax, out=hs)
        array = np.round((hs - offset) / scale).astype('uint8')
        array[hs.mask] = 255
        colortable = stfmt.format_colortable('matplotlib_jet', vmin=vmin, vmax=vmax,
                                             vmin_pal=vmin_pal, vmax_pal=vmax_pal)
        band.append({'array': array,
                     'scale': scale,
                     'offset': offset,
                     'description': metadata['parameter'],
                     'unittype': 'm',
                     'nodatavalue': 255,
                     'parameter_range': [vmin, vmax],
                     'colortable': colortable})
        # Write geotiff
        print 'Write geotiff'
        tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
        stfmt.write_geotiff(tifffile, metadata, geolocation, band)

        ### Partitions ###
        phslst = [ww3['phs'+str(i)][itime, :, :] for i in range(ww3['npart'])]
        phs = np.ma.dstack(phslst)
        ptplst = [ww3['ptp'+str(i)][itime, :, :] for i in range(ww3['npart'])]
        ptp = np.ma.dstack(ptplst)
        pdirlst = [ww3['pdir'+str(i)][itime, :, :] for i in range(ww3['npart'])]
        pdir = np.ma.dstack(pdirlst)
        # Reorder partitions by HS -> Keep WW3 order
        # phs.data[phs.mask] = -1000 # make sure masked values don't interfere with sorting
        # index = np.ogrid[:phs.shape[0], :phs.shape[1], :phs.shape[2]]
        # index[2] = (-phs).argsort(axis=2, kind='mergesort')
        # phs = phs[index]
        # ptp = ptp[index]
        # pdir = pdir[index]
        # pdir from_direction -> to_direction
        pdir = np.mod(pdir + 180.0, 360.0)
        # pdir clockwise from north -> counter clockwise from east
        pdir = np.mod(90.0 - pdir, 360.0)
        # Write each partition in a geotiff
        for i in range(ww3['npart']):
            # Update metadata
            lpartnum = 'partition '+str(i)
            spartnum = 'part'+str(i)
            metadata['product_name'] = 'WW3_model_wave_'+spartnum
            if v2 == True:
                metadata['product_name'] += '_v2'
            metadata['name'] = basename + '-' + spartnum
            metadata['parameter'] = ['wave significant height '+lpartnum,
                                     'wave peak period '+lpartnum,
                                     'wave mean direction '+lpartnum]
            # Make bands
            band = []
            iphs, iptp, ipdir = phs[:, :, i], ptp[:, :, i], pdir[:, :, i]
            # HS
            #_vmin, _vmax = 0.0, 25.4
            offset, scale = vmin, (vmax - vmin) / 254.0
            np.clip(iphs, vmin, vmax, out=iphs)
            array = np.round((iphs - offset) / scale).astype('uint8')
            array[iphs.mask] = 255
            band.append({'array': array,
                         'scale': scale,
                         'offset': offset,
                         'description': metadata['parameter'][0],
                         'unittype': 'm',
                         'nodatavalue': 255,
                         'parameter_range': [vmin, vmax]})
            # Period
            _vmin, _vmax = 0.0, 25.4
            offset, scale = _vmin, (_vmax - _vmin) / 254.0
            np.clip(iptp, _vmin, _vmax, out=iptp)
            array = np.round((iptp - offset) / scale).astype('uint8')
            array[iptp.mask] = 255
            band.append({'array': array,
                         'scale': scale,
                         'offset': offset,
                         'description': metadata['parameter'][1],
                         'unittype': 's',
                         'nodatavalue': 255,
                         'parameter_range': [_vmin, _vmax]})
            # Direction
            _vmin, _vmax = 0.0, 360.0
            offset, scale = _vmin, (_vmax - _vmin) / 254.0
            np.clip(ipdir, _vmin, _vmax, out=ipdir)
            array = np.round((ipdir - offset) / scale).astype('uint8')
            array[ipdir.mask] = 255
            band.append({'array': array,
                         'scale': scale,
                         'offset': offset,
                         'description': metadata['parameter'][2],
                         'unittype': 'degree',
                         'nodatavalue': 255,
                         'parameter_range': [vmin, vmax]})
            # Write geotiff
            print 'Write geotiff'
            tifffile = stfmt.format_tifffilename(outdir, metadata, create_dir=True)
            stfmt.write_geotiff(tifffile, metadata, geolocation, band,
                                drv_opts=['PHOTOMETRIC=MINISBLACK'])
