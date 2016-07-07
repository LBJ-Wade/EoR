# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
from astropy.io import fits
from astropy import units
from os.path import join
import matplotlib.pyplot as plt
from scipy import optimize
from astropy import constants as const


if __name__ == '__main__':
    # http://goo.gl/Yh2RwR
    # http://goo.gl/ra2eWq

    beam_sigma = 50 * units.arcsec
    beam_area = 2 * np.pi * beam_sigma**2
    freq = 100 * units.MHz
    equiv = units.brightness_temperature(beam_area, freq)
    units.Jy.to(units.K, equivalencies=equiv)
    print('1jy -> K = ', (1 * units.Jy).to(units.K, equivalencies=equiv))

    conv = const.c.value**2 / (2 * const.k_B.value * freq.to(units.Hz).value**2)
    print('beam area', beam_area, beam_area.to(units.sr))
    print('conv', conv)

    # 1 Jy = 1e-26 W (m^-2) (Hz^-2) -> K
    print((1.e-26 / beam_area.to(units.sr).value) * conv)

    # TODO(BM)
    # 1 load cube
    # 2 load beam areas
    # 3 convert to K and save out cube. update BUNITS header.

    # psf_cube = fits.getdata(join('results',
    #                              'noise_050.0-070.0MHz_100kHz_5h_60s_natural',
    #                              'noise_l_cut_10_230_natural.fits'))
    #
    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(111)
    # im = ax.imshow(psf_cube[0, :, :], interpolation='nearest', origin='lower',
    #                cmap='inferno')
    # ax.figure.colorbar(im, ax=ax)
    # plt.show()

