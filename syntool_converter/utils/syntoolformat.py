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

import gdal
import osr
import os
import numpy as np
from datetime import datetime, timedelta
from syntool_converter.utils.gmtColormap import gmtColormap
from netCDF4 import Dataset
import pyproj


TIMEFMT = '%Y-%m-%dT%H:%M:%S'
NP2GDAL = {'uint8':gdal.GDT_Byte, 'uint16':gdal.GDT_UInt16,
           'int16':gdal.GDT_Int16, 'uint32':gdal.GDT_UInt32,
           'int32':gdal.GDT_Int32, 'float32':gdal.GDT_Float32,
           'float64':gdal.GDT_Float64}
COL2GDAL = {'gray': gdal.GCI_GrayIndex, 'palette': gdal.GCI_PaletteIndex,
            'alpha': gdal.GCI_AlphaBand, 'red': gdal.GCI_RedBand,
            'green': gdal.GCI_GreenBand, 'blue': gdal.GCI_BlueBand,
            'undefined': gdal.GCI_Undefined}
CMAPS = [ 'cerbere_wind'
        , 'cerbere_medspiration'
        , 'rainbow'
        , 'noaa_wind'
        , 'sea_ice_conc'
        , 'matplotlib_jet'
        , 'matplotlib_jet_r'
        , 'doppler'
        , 'gmt_panoply'
        , 'yellow_red'
        , 'pesket'
        , 'gebco'
        , 'ibcso'
        , 'chla_jet'
        ]

def format_time_and_range(start_time, stop_time, units='ms'):
    """
    """
    # Compute middle time, round to second and format
    middle_time = start_time + (stop_time-start_time)/2
    microsecond = middle_time.microsecond
    middle_time = middle_time.replace(microsecond=0)
    if microsecond > 500000:
        middle_time += timedelta(seconds=1)
    dtime = format_time(middle_time)
    # Compute delta times and format
    delta_start = (start_time-middle_time).total_seconds()
    delta_stop = (stop_time-middle_time).total_seconds()
    fac = {'h':1./3600, 'm':1./60, 's':1, 'ms':1000}[units]
    fmt = '%+i'+units
    time_range = [fmt % np.round(delta_start*fac),
                  fmt % np.round(delta_stop*fac)]
    # Return time and range
    return (dtime, time_range)


def format_time(time):
    """
    """
    return time.strftime(TIMEFMT)


def format_tifffilename(outdir, metadata, create_dir=False):
    """
    """
    if os.environ.get('SYNTOOL_CONVERTER_KEEP_OUTPUT', '0') == '1':
        tiffdir = outdir
    else:
        pdname = metadata['product_name'].lower().replace(' ', '_')
        tiffdir = os.path.join(outdir, pdname)
        if 'datagroup' in metadata:
            tiffdir = os.path.join(tiffdir, metadata['datagroup'])
    if create_dir == True and not os.path.exists(tiffdir):
        try:
            os.makedirs(tiffdir)
        except OSError:
            # Check that the directory was not created by another process in
            # the meantime.
            if not os.path.isdir(tiffdir):
                raise OSError('Failed to create {}'.format(tiffdir))
    dsname = metadata['name'].replace(' ', '_')
    tifffile = os.path.join(tiffdir, dsname+'.tiff')
    return tifffile


def format_gdalgcps(x, y, z, pix, lin):
    """
    """
    sizes = [x.size, y.size, z.size, pix.size, lin.size]
    if min(sizes) != max(sizes):
        raise Exception('All inputs must be of same size.')
    gcps = []
    for i, j, k, l, m in zip(x.flat, y.flat, z.flat, pix.flat, lin.flat):
        gcps.append(gdal.GCP(float(i), float(j), float(k), float(l), float(m)))
    return gcps


def format_gdalprojection(geogcs='WGS84'):
    """
    """
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS(geogcs)
    return srs.ExportToWkt()


def write_geotiff(tifffile, metadata, geolocation, band, drv_opts=[]):
    """
    """
    ### driver options examples
    #drv_opts = []
    #drv_opts = ['PHOTOMETRIC=MINISBLACK']
    #drv_opts = ['PHOTOMETRIC=PALETTE']
    #drv_opts = ['RGB']
    # Create dataset
    nband = len(band)
    shape = band[0]['array'].shape
    datatype = NP2GDAL[str(band[0]['array'].dtype)]
    drv = gdal.GetDriverByName('GTiff')
    dset = drv.Create(tifffile, shape[1], shape[0], nband, datatype, drv_opts)
    # Write metadata
    for name, value in metadata.iteritems():
        if isinstance(value, list):
            strvalue = []
            for listvalue in iter(value):
                if isinstance(listvalue, basestring):
                    strvalue.append('"'+str(listvalue)+'"')
                else:
                    strvalue.append(str(listvalue))
            strvalue = ' '.join(strvalue)
        elif isinstance(value, basestring):
            strvalue = '"'+str(value)+'"'
        else:
            strvalue = str(value)
        dset.SetMetadataItem(name, strvalue)
    # Write geolocation
    if 'geotransform' in geolocation.keys():
        dset.SetProjection(geolocation['projection'])
        dset.SetGeoTransform(geolocation['geotransform'])
    elif 'gcps' in geolocation.keys():
        dset.SetGCPs(geolocation['gcps'], geolocation['projection'])
    else:
        raise Exception('Need geotransform or gcps')
    # Write band(s)
    for ibnd in range(nband):
        bnd = dset.GetRasterBand(int(ibnd+1))
        for name, value in band[ibnd].iteritems():
            if name == 'array':
                bnd.WriteArray(value)
            elif name == 'scale':
                bnd.SetScale(value)
            elif name == 'offset':
                bnd.SetOffset(value)
            elif name == 'description':
                bnd.SetDescription('"'+value+'"')
            elif name == 'unittype':
                bnd.SetUnitType('"'+value+'"')
            elif name == 'nodatavalue':
                bnd.SetNoDataValue(value)
            elif name == 'valid_range': # Use parameter_range instead !
                bnd.SetMetadataItem(name, str(value[0])+' '+str(value[1]))
            elif name == 'parameter_range':
                bnd.SetMetadataItem(name, str(value[0])+' '+str(value[1]))
            elif name == 'colortable':
                bnd.SetColorTable(value)
            elif name == 'colorinterp':
                bnd.SetColorInterpretation(COL2GDAL[value])
    # Close dataset
    dset = None


def get_shape(dataset, npts=10):
    """
    """
    transformer = gdal.Transformer(dataset, None, ['MAX_GCP_ORDER=-1'])
    npix = dataset.RasterXSize
    nlin = dataset.RasterYSize
    pix = np.hstack((np.round(np.linspace(0, npix, num=npts)),
                     np.repeat(npix, npts),
                     np.round(np.linspace(npix, 0, num=npts)),
                     np.repeat(0, npts)))
    lin = np.hstack((np.repeat(0, npts),
                     np.round(np.linspace(0, nlin, num=npts)),
                     np.repeat(nlin, npts),
                     np.round(np.linspace(nlin, 0, num=npts))))
    pixlin = np.vstack((pix, lin)).transpose()
    lonlat = np.array(transformer.TransformPoints(0, pixlin)[0])
    return (lonlat[:, 0], lonlat[:, 1])


# def _test(tifffile):
#     """
#     """
#     dataset = gdal.Open(tifffile, gdal.GA_ReadOnly)
#     (lon, lat) = get_shape(dataset)
#     srs_src = osr.SpatialReference()
#     #srs_src.SetWellKnownGeogCS('WGS84')
#     srs_src.ImportFromEPSG(4326)
#     srs_dst = osr.SpatialReference()
#     srs_dst.ImportFromEPSG(900913)
#     transformer = osr.CoordinateTransformation(srs_src, srs_dst)
#     lonlat = np.vstack((lon, lat)).transpose()
#     xy = transformer.TransformPoints(lonlat)
#     pdb.set_trace()


def write_pngkml_proj(tifffile):
    """
    """
    # Get infos from tiff
    basefile = os.path.splitext(tifffile)[0]
    dset = gdal.Open(tifffile, gdal.GA_ReadOnly)
    product_name = dset.GetMetadataItem('product_name')[1:-1]
    timestamp = dset.GetMetadataItem('datetime')[1:-1]
    (shplon, shplat) = get_shape(dset)
    # bnd1 = dset.GetRasterBand(1)
    # nodatavalue = bnd1.GetNoDataValue()
    dset = None
    # Convert tiff to projected tiff
    projtifffile = basefile+'-proj.tiff'
    ### googlemap EPSG:3857 (mercator)
    ### googleearth EPSG:4326 (equirectangular)
    ### syntool EPSG:900913 (Google Maps Global Mercator)
    # cmd = 'gdalwarp -t_srs EPSG:4326 -order 3 -r bilinear '+ \
    #       '-dstnodata '+str(nodatavalue)+' '+ \
    #       '-overwrite '+tifffile+' '+projtifffile
    ##'-tr '+str(spacing[1]*2)+' '+str(spacing[0]*2)+' '+ \
    cmd = 'gdalwarp -t_srs EPSG:4326 -tps -r bilinear -dstalpha '+ \
          '-overwrite '+tifffile+' '+projtifffile
    os.system(cmd)
    # Convert projected tiff to projected png
    projpngfile = basefile+'-proj.png'
    # cmd = 'gdal_translate -a_nodata '+str(nodatavalue)+' -of PNG '+ \
    #       projtifffile+' '+projpngfile
    cmd = 'gdal_translate -of PNG '+projtifffile+' '+projpngfile
    os.system(cmd)
    os.remove(projtifffile)
    os.remove(projpngfile+'.aux.xml')
    # Create kml
    projkmlfile = basefile+'-proj.kml'
    kml1 = ['<?xml version="1.0" encoding="UTF-8"?>',
            '<kml xmlns="http://earth.google.com/kml/2.0" xmlns:gx="http://www.google.com/kml/ext/2.2">',
            '  <Document>',
            '    <name>'+product_name+'</name>',
            '    <description>'+timestamp+'</description>',
            '    <GroundOverlay>',
            '      <name>Image layer</name>',
            '      <TimeStamp>',
            '        <when>'+timestamp+'</when>',
            '      </TimeStamp>',
            '      <Icon>',
            '        <href>'+os.path.basename(projpngfile)+'</href>',
            '      </Icon>',
            '      <LatLonBox>',
            '        <north>'+str(shplat.max())+'</north>',
            '        <south>'+str(shplat.min())+'</south>',
            '        <west>'+str(shplon.min())+'</west>',
            '        <east>'+str(shplon.max())+'</east>',
            '      </LatLonBox>',
            '    </GroundOverlay>',
            '    <Placemark>',
            '      <name>Bounding polygon</name>',
            '      <styleUrl>#m_ylw-pushpin</styleUrl>',
            '      <Polygon>',
            '        <tessellate>1</tessellate>',
            '        <outerBoundaryIs>',
            '          <LinearRing>',
            '            <coordinates>']
    kml2 = []
    for i in np.arange(shplon.size):
        kml2.append('              '+str(shplon[i])+','+str(shplat[i])+',0')
    kml3 = ['            </coordinates>',
            '          </LinearRing>',
            '        </outerBoundaryIs>',
            '      </Polygon>',
            '    </Placemark>',
            '  </Document>',
            '</kml>']
    kmllines = kml1+kml2+kml3
    fil = open(projkmlfile, "w")
    for kmlline in kmllines:
        fil.write(kmlline+'\n')
    fil.close()


def load_colormap(name):
    """
    """
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    import matplotlib.cm
    palettes = { 'cerbere_wind': 'wind.pal'
               , 'cerbere_medspiration': 'medspiration.rgb'
               , 'rainbow': 'rainbow.rgb'
               , 'gmt_panoply': 'GMT_panoply.cpt' 
               , 'gebco': 'GMT_gebco.cpt'
               ,'ibcso': 'ibcso-bath.cpt'
               ,'chla_jet': 'chla_jet.cpt'
               }
    fname = palettes.get(name, None)
    if fname is not None:
        if os.path.splitext(fname)[1] == '.cpt':
            cpt_file = os.path.join(os.path.dirname(__file__), 'palette', fname)
            colordict = gmtColormap(cpt_file)
            colormap = LinearSegmentedColormap('custom', colordict)
            nodata_color = (0, 0, 0)
        else:
            rgb_file = os.path.join(os.path.dirname(__file__), 'palette', fname)
            colors = []
            for line in open(rgb_file):
                colors.append([float(col)/255 for col in line.split()])
            colormap = ListedColormap(colors, name='custom')
            #colormap = LinearSegmentedColormap.from_list('custom', colors, N=256)
            nodata_color = (0, 0, 0)
    elif name == 'noaa_wind':
        ### From cpt file :
        # COLOR_MODEL = RGB
        # 0   0   0   0	  0   0   0   0
        # 1 127 127 127	  5  44  44  44
        # 5   0 255 255	 15   0   0 255
        # 15   0 255   0	 20   0 115   0
        # 20 255 255  43	 25 255 186   0
        # 25 255 145  12	 30 232  33   0
        # 30 255  43  70	 35 255   0   0
        # 35 205 127  78	 40  71  12  12
        # 40 255   0 255	 45  98   0 154
        # 50  50   0 100	 50  50   0 100
        # B	80	0	100
        # F	255	0	0
        # N	255	255	255
        from_cpt = [(0, 0, 0, 0), (0, 0, 0, 0),
                    (1, 127, 127, 127), (1, 127, 127, 127),
                    (5, 44, 44, 44), (5, 0, 255, 255),
                    (15, 0, 0, 255), (15, 0, 255, 0),
                    (20, 0, 115, 0), (20, 255, 255, 43),
                    (25, 255, 186, 0), (25, 255, 145, 12),
                    (30, 232, 33, 0), (30, 255, 43, 70),
                    (35, 255, 0, 0), (35, 205, 127, 78),
                    (40, 71, 12, 12), (40, 255, 0, 255),
                    (45, 98, 0, 154), (45, 98, 0, 154),
                    (50, 50, 0, 100), (50, 50, 0, 100)]
        cdict = {'red':[], 'green':[], 'blue':[]}
        for i in range(0, len(from_cpt), 2):
            n1, r1, g1, b1 = from_cpt[i]
            n2, r2, g2, b2 = from_cpt[i+1]
            cdict['red'].append((n1/50., r1/255., r2/255.))
            cdict['green'].append((n1/50., g1/255., g2/255.))
            cdict['blue'].append((n1/50., b1/255., b2/255.))
        colormap = LinearSegmentedColormap('custom', cdict, N=256)
        nodata_color = (255, 255, 255)
    elif name == 'sea_ice_conc':
        colors = ["#08316C", "#ffffff"] # RGB : [(8, 49, 108), (255, 255, 255)]
        colormap = LinearSegmentedColormap.from_list('custom', colors, N=256)
        nodata_color = (0, 0, 0)
    elif name == 'matplotlib_jet':
        colormap = matplotlib.cm.jet
        nodata_color = (255, 255, 255)
    elif name == 'matplotlib_jet_r':
        colormap = matplotlib.cm.jet_r
        nodata_color = (255, 255, 255)
    elif name == 'doppler':
        xval = [-2.5, -2.0, -1.6, -0.8, -0.3, 0, 0.3, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5]
        red = [255, 220, 128, 12, 50, 255, 90, 255, 255, 255, 255, 255, 160]
        green = [0, 0, 30, 120, 200, 255, 255, 255, 255, 170, 85, 0, 0]
        blue = [128, 220, 200, 238, 235, 255, 30, 0, 0, 0, 0, 0, 0]
        segmentdata = {'red':[], 'green':[], 'blue':[]}
        for x, r, g, b in zip(xval, red, green, blue):
            xnorm = (x+2.5)/5
            rnorm = r/255.
            gnorm = g/255.
            bnorm = b/255.
            segmentdata['red'].append((xnorm, rnorm, rnorm))
            segmentdata['green'].append((xnorm, gnorm, gnorm))
            segmentdata['blue'].append((xnorm, bnorm, bnorm))
        colormap = LinearSegmentedColormap('custom', segmentdata)
        nodata_color = (0, 0, 0)
    elif name == 'yellow_red':
        colors = ["#FFFF00", "#FF0000"]
        colormap = LinearSegmentedColormap.from_list('custom', colors)
        nodata_color = (0, 0, 0)
    elif name == 'pesket':
        # edges = [-1, 6, 13, 19, 26, 32, 39, 45, 52, 59, 65, 72, 78, 85, 91, 98,
        #          105, 111, 118, 124, 131, 137, 144, 150, 157, 164, 170, 177,
        #          183, 190, 196, 203, 210, 216, 223, 229, 236, 242, 249, 254]
        edges = [0., 6.5, 13.5, 19.5, 26.5, 32.5, 39.5, 45.5, 52.5, 59.5,
                 65.5, 72.5, 78.5, 85.5, 91.5, 98.5, 105.5, 111.5, 118.5,
                 124.5, 131.5, 137.5, 144.5, 150.5, 157.5, 164.5, 170.5,
                 177.5, 183.5, 190.5, 196.5, 203.5, 210.5, 216.5, 223.5,
                 229.5, 236.5, 242.5, 249.5, 254.]
        enrm = [float(e) / 254 for e in edges]
        rnrm = [0.00, 0.00, 0.43, 0.53, 0.67, 0.50, 0.60, 0.72, 0.60, 0.43,
                0.12, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.78,
                0.77, 0.87, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                1.00, 1.00, 1.00, 1.00, 0.87, 0.75, 0.60, 0.53, 0.40]
        gnrm = [0.00, 0.00, 0.45, 0.55, 0.70, 0.80, 0.87, 0.95, 0.97, 0.97,
                0.97, 0.97, 0.93, 0.55, 0.60, 0.75, 0.88, 1.00, 0.97, 0.96,
                0.88, 0.95, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.83, 0.75,
                0.65, 0.58, 0.50, 0.30, 0.00, 0.00, 0.00, 0.38, 0.40]
        bnrm = [0.55, 0.90, 0.82, 0.92, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
                0.90, 0.80, 0.82, 0.53, 0.00, 0.00, 0.00, 0.00, 0.48, 0.60,
                0.40, 0.48, 0.80, 0.62, 0.42, 0.00, 0.72, 0.65, 0.62, 0.50,
                0.38, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.40]
        segmentdata = {'red':[], 'green':[], 'blue':[]}
        for e1, e2, r, g, b in zip(enrm[0:-1], enrm[1:], rnrm, gnrm, bnrm):
            segmentdata['red'].append((e1, r, r))
            segmentdata['red'].append((e2, r, r))
            segmentdata['green'].append((e1, g, g))
            segmentdata['green'].append((e2, g, g))
            segmentdata['blue'].append((e1, b, b))
            segmentdata['blue'].append((e2, b, b))
        colormap = LinearSegmentedColormap('custom', segmentdata)
        nodata_color = (255, 255, 255)
    else:
        raise Exception('Unknown colormap name.')
    return [colormap, nodata_color]


def format_colortable(name, vmin=0., vmax=1., vmin_pal=0., vmax_pal=1.,
                      index_min=0, index_max=254, index_nodata=255):
    """
    """
    colormap, nodata_color = load_colormap(name)
    norm_min = (vmin - vmin_pal) / float((vmax_pal - vmin_pal))
    norm_max = (vmax - vmin_pal) / float((vmax_pal - vmin_pal))
    ncols = int(index_max) - int(index_min) + 1
    colors = colormap(np.linspace(norm_min, norm_max, num=ncols))
    colors = np.round(colors*255)
    entries = [(int(c[0]), int(c[1]), int(c[2])) for c in colors]
    colortable = gdal.ColorTable()
    for index in range(int(index_min), int(index_max)+1):
        colortable.SetColorEntry(index, entries[index-int(index_min)])
    if index_nodata != None:
        colortable.SetColorEntry(index_nodata, nodata_color)
    return colortable


def format_ncfilename(outdir, metadata, create_dir=False):
    """
    """
    pdname = metadata['product_name'].lower().replace(' ', '_')
    dsname = metadata['name'].replace(' ', '_')
    if os.environ.get('SYNTOOL_CONVERTER_KEEP_OUTPUT', '0') == '1':
        ncdir = outdir
    else:
        ncdir = os.path.join(outdir, pdname, dsname)
        # if 'datagroup' in metadata:
        #     ncdir = os.path.join(ncdir, metadata['datagroup'])
    if create_dir == True and not os.path.exists(ncdir):
        try:
            os.makedirs(ncdir)
        except OSError:
            # Check that the directory was not created by another process in
            # the meantime.
            if not os.path.isdir(ncdir):
                raise OSError('Failed to create {}'.format(ncdir))
    ncfile = os.path.join(ncdir, '{}_idf_00.nc'.format(dsname))
    return ncfile


def write_netcdf(ncfile, metadata, geolocation, band, model,
                 ngcps=None, dgcps=None):
    """
    """
    nband = len(band)
    nlin, npix = band[0]['array'].shape

    # GCPs
    if ngcps is None and dgcps is None:
        raise Exception('Need number of GCPs or pixel distance between GCPs.')
    elif ngcps is None:
        ngcps = (np.ceil(nlin / float(dgcps[0])).astype('int') + 1,
                 np.ceil(npix / float(dgcps[1])).astype('int') + 1)
    elif dgcps is None:
        dgcps = (np.ceil(nlin / (ngcps[0] - 1.)).astype('int'),
                 np.ceil(npix / (ngcps[1] - 1.)).astype('int'))
    gcplin = np.concatenate((np.arange(ngcps[0] - 1) * dgcps[0], [nlin]))
    gcppix = np.concatenate((np.arange(ngcps[1] - 1) * dgcps[1], [npix]))
    if 'geotransform' in geolocation.keys():
        x0, dxx, dxy, y0, dyx, dyy = geolocation['geotransform']
        gcpy = y0 + dyx * gcppix[np.newaxis, :] + dyy * gcplin[:, np.newaxis]
        gcpx = x0 + dxx * gcppix[np.newaxis, :] + dxy * gcplin[:, np.newaxis]
        srs = osr.SpatialReference()
        srs.ImportFromWkt(geolocation['projection'])
        srs4326 = osr.SpatialReference()
        srs4326.ImportFromEPSG(4326)
        if srs.IsSame(srs4326) == 1:
            gcplat = gcpy
            gcplon = gcpx
        else:
            proj = pyproj.Proj(srs.ExportToProj4())
            gcplon, gcplat = proj(gcpx, gcpy, inverse=True)
    elif 'gcps' in geolocation.keys():
        drv = gdal.GetDriverByName('MEM')
        dset = drv.Create('tmp', npix, nlin)
        dset.SetGCPs(geolocation['gcps'], geolocation['projection'])
        options = ['MAX_GCP_ORDER=-1']
        transformer = gdal.Transformer(dset, None, options)
        _gcplin = np.tile(gcplin[:, np.newaxis], (1, ngcps[1]))
        _gcppix = np.tile(gcppix[np.newaxis, :], (ngcps[0], 1))
        inputs = np.vstack((_gcppix.flatten(), _gcplin.flatten())).transpose()
        outputs = np.array(transformer.TransformPoints(0, inputs)[0])[:, 0:2]
        gcplat = outputs[:, 1].reshape(ngcps)
        gcplon = outputs[:, 0].reshape(ngcps)
        dset = None

    # Time
    time = datetime.strptime(metadata['datetime'], TIMEFMT)
    time_range = []
    for dts in metadata['time_range']:
        if dts.endswith('ms'):
            timed = timedelta(milliseconds=int(dts[:-2]))
        elif dts.endswith('s'):
            timed = timedelta(seconds=int(dts[:-1]))
        elif dts.endswith('m'):
            timed = timedelta(minutes=int(dts[:-1]))
        elif dts.endswith('h'):
            timed = timedelta(hours=int(dts[:-1]))
        elif dts.endswith('d'):
            timed = timedelta(days=int(dts[:-1]))
        else:
            raise Exception()
        time_range.append(time + timed)

    # Model
    if model == 'grid_lonlat':
        dim1name = 'lat'
        dim2name = 'lon'
        lonlat1d = True
        unlimtime = True
    elif model == 'grid_proj':
        dim1name = 'y'
        dim2name = 'x'
        lonlat1d = False
        unlimtime = True
    elif model == 'swath':
        dim1name = 'row'
        dim2name = 'cell'
        lonlat1d = False
        unlimtime = False
    else:
        raise Exception('')
    dim1gcpname = '{}_gcp'.format(dim1name)
    dim2gcpname = '{}_gcp'.format(dim2name)

    # Write
    dset = Dataset(ncfile, mode='w', format='NETCDF4', clobber=True)
    ## Dimensions
    if unlimtime:
        _dimtime = dset.createDimension('time', None)
    else:
        _dimtime = dset.createDimension('time', 1)
    _dim1 = dset.createDimension(dim1name, nlin)
    _dim2 = dset.createDimension(dim2name, npix)
    _dim1gcp = dset.createDimension(dim1gcpname, ngcps[0])
    _dim2gcp = dset.createDimension(dim2gcpname, ngcps[1])
    ## Variables
    _time = dset.createVariable('time', 'f8', ('time', ))
    _time.long_name = 'time'
    _time.standard_name = 'time'
    _time.units = 'seconds since 1970-01-01T00:00:00.000000Z'
    _time.calendar = 'standard'
    _time[:] = (time - datetime(1970, 1, 1)).total_seconds()
    if lonlat1d:
        _latgcp = dset.createVariable('lat_gcp', 'f4', (dim1gcpname,))
        _latgcp[:] = gcplat[:, 0].astype('float32')
    else:
        _latgcp = dset.createVariable('lat_gcp', 'f4', (dim1gcpname, dim2gcpname))
        _latgcp[:] = gcplat.astype('float32')
    _latgcp.long_name = 'ground control points latitude'
    _latgcp.standard_name = 'latitude'
    _latgcp.units = 'degrees_north'
    if lonlat1d:
        _longcp = dset.createVariable('lon_gcp', 'f4', (dim2gcpname,))
        _longcp[:] = gcplon[0, :].astype('float32')
    else:
        _longcp = dset.createVariable('lon_gcp', 'f4', (dim1gcpname, dim2gcpname))
        _longcp[:] = gcplon.astype('float32')
    _longcp.long_name = 'ground control points longitude'
    _longcp.standard_name = 'longitude'
    _longcp.units = 'degrees_east'
    _indexdim1gcp = dset.createVariable('index_{}_gcp'.format(dim1name), 'i4', (dim1gcpname,))
    _indexdim1gcp[:] = gcplin.astype('int32')
    _indexdim1gcp.long_name = 'index of ground control points in {} dimension'.format(dim1name)
    _indexdim1gcp.comment = 'index goes from 0 (start of first pixel) to dimension value '\
                            '(end of last pixel)'
    _indexdim2gcp = dset.createVariable('index_{}_gcp'.format(dim2name), 'i4', (dim2gcpname,))
    _indexdim2gcp[:] = gcppix.astype('int32')
    _indexdim2gcp.long_name = 'index of ground control points in {} dimension'.format(dim2name)
    _indexdim2gcp.comment = 'index goes from 0 (start of first pixel) to dimension value '\
                            '(end of last pixel)'
    for ibnd in range(nband):
        if 'name' in band[ibnd]:
            varname = band[ibnd]['name']
        elif 'description' in band[ibnd]:
            varname = band[ibnd]['description'].replace(' ', '_')
        else:
            varname = 'value_{}'.format(ibnd)
        dtype = band[ibnd]['array'].dtype
        fill_value = None
        if 'nodatavalue' in band[ibnd]:
            fill_value = dtype.type(band[ibnd]['nodatavalue'])
        _value = dset.createVariable(varname, dtype, ('time', dim1name, dim2name),
                                     fill_value=fill_value)
        _value[:] = band[ibnd]['array'][np.newaxis, :, :]
        if 'long_name' in band[ibnd]:
            _value.long_name = band[ibnd]['long_name']
        if 'standard_name' in band[ibnd]:
            _value.standard_name = band[ibnd]['standard_name']
        if 'comment' in band[ibnd]:
            _value.comment = band[ibnd]['comment']
        # if 'description' in band[ibnd]:
        #     _value.long_name = band[ibnd]['description']
        if 'unittype' in band[ibnd]:
            _value.units = band[ibnd]['unittype']
        if 'offset' in band[ibnd]:
            _value.add_offset = np.float32(band[ibnd]['offset'])
        if 'scale' in band[ibnd]:
            _value.scale_factor = np.float32(band[ibnd]['scale'])
        if 'valid_range' in band[ibnd]:
            _value.valid_min = dtype.type(band[ibnd]['valid_range'][0])
            _value.valid_max = dtype.type(band[ibnd]['valid_range'][1])
        if 'parameter_range' in band[ibnd]:
            offset = 0.
            if 'offset' in band[ibnd]:
                offset = band[ibnd]['offset']
            scale = 1.
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
    dset.time_coverage_start = time_range[0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    dset.time_coverage_end = time_range[1].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    dset.idf_subsampling_factor = np.int32(0)
    if 'spatial_resolution' in metadata:
        dset.idf_spatial_resolution = np.float32(metadata['spatial_resolution'])
    else:
        raise Exception('IDF spatial resolution ?')
    dset.idf_spatial_resolution_units = 'm'
    dset.close()

