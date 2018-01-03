# coding: utf-8

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
import logging
import numpy as np

from matplotlib import pyplot
import matplotlib as mpl

import syntoolformat as stfmt

def identity(value):
    return value

def toKnots(value):
    return value * 1.94384

def toCelsius(value):
    return value - 273.15

logger = logging.getLogger(__name__)

CONVERT_METHODS = {
    'identity': identity,
    'toKnots': toKnots,
    'toCelsius': toCelsius
}

def create_colorbar_image(output_dir, product_config):
    fig = pyplot.figure(figsize=(4, 0.7), dpi=80)
    axis = fig.add_axes([0.05, 0.625, 0.9, 0.30])
    fig.patch.set_facecolor('white')

    cmap, _ = stfmt.load_colormap(product_config['colorbar'])
    convert = CONVERT_METHODS[product_config['convert']]
    relative = product_config.get('relative', False)
    discrete = product_config.get('discrete', False)
    if not relative:
        vmin_pal = convert(product_config['min'])
        vmax_pal = convert(product_config['max'])
    else:
        vmin_pal = 0.
        vmax_pal = 1.
    norm = mpl.colors.Normalize(vmin=vmin_pal, vmax=vmax_pal)
    nbins = product_config.get('nbins', None)

    if not discrete:
        cb = mpl.colorbar.ColorbarBase(axis, cmap=cmap, norm=norm,
                                   orientation='horizontal')
    else:
        bounds = product_config.get('bounds',
                                    np.linspace(vmin_pal, vmax_pal, num=nbins))
        cb = mpl.colorbar.ColorbarBase(axis, cmap=cmap, norm=norm,
                                       boundaries=[0]+ bounds + [0],
                                       extend='both', orientation='horizontal',
                                       spacing='proportional')

    if nbins is not None:
        locator = mpl.ticker.MaxNLocator(nbins=nbins)
        cb.locator = locator
        cb.update_ticks()

    if relative:
        nticklabels = len(cb.ax.get_xticklabels())
        ticklabels = ['min'] + [''] * (nticklabels - 2) + ['max']
        cb.ax.set_xticklabels(ticklabels)

    label = product_config["variable"]
    if product_config["units"]:
        units = product_config['units']
        label = u'{} ({})'.format(label, units)
    cb.set_label(label, labelpad=-0.25)

    file_name = u'{}_colorbar.png'.format(product_config['product'])
    image_path = os.path.join( output_dir, file_name)
    logger.info('Create colorbar: ' + image_path)
    pyplot.savefig(image_path, dpi=80, facecolor=fig.get_facecolor(), edgecolor='none')
    pyplot.close(fig)

