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
from scipy.ndimage.filters import convolve
from scipy.signal.windows import gaussian


def denoise_sst(sst, bt11, kernel='boxcar', width=20, show=False):
    """
    """
    if kernel == 'boxcar':
        ker = np.ones((width, width))
    elif kernel == 'gaussian':
        gau = gaussian(3*width, width/2.)
        ker = np.outer(gau, gau)
    else:
        raise Exception('Unknown kernel')
    mask = ma.getmaskarray(sst) | ma.getmaskarray(bt11)
    sst.data[mask] = 0.
    bt11.data[mask] = 0.
    corr = convolve((sst - 273.15) * 0.9 + 273.15 - bt11, ker)
    #corr = convolve(sst - bt11, ker)
    #corr = convolve((sst - 273) * 0.9 + 273 - bt11, ker)
    norm = convolve(1. - mask, ker)
    corr[mask] = 0
    corr[~mask] /= norm[~mask]
    denoised_sst = bt11 + corr
    denoised_sst.mask = mask
    if show == True:
        import matplotlib.pyplot as plt
        vmin = min([bt11.min(), sst.min(), denoised_sst.min()])
        vmax = max([bt11.max(), sst.max(), denoised_sst.max()])
        corr = ma.MaskedArray(corr, mask=mask)
        for ivar, var in enumerate([sst, bt11, corr, denoised_sst]):
            plt.figure()
            if ivar == 2:
                plt.imshow(var, interpolation='nearest')
            else:
                plt.imshow(var, vmin=vmin, vmax=vmax, interpolation='nearest')
            plt.colorbar()
        plt.show()
    return denoised_sst
