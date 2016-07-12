# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
from astropy.io import fits
from astropy import units
from os.path import join
import matplotlib.pyplot as plt
from scipy import optimize
from math import radians, pi
from astropy import constants as const
import os


def main():
    for freq0 in (50e6, 120e6, 200e6):
        for weight in ('uniform', 'natural'):
            print(freq0/1e6,  weight)

            freq1 = freq0 + 20e6
            fit_file = join('results', 'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
                            (freq0/1e6, freq1/1e6, weight),
                            'psf_fit', 'fit_data.txt')
            if freq0 == 50e6:
                b_min_lambda = 10
                b_max_lambda = 230
            elif freq0 == 120e6:
                b_min_lambda = 20
                b_max_lambda = 550
            elif freq0 == 200e6:
                b_min_lambda = 30
                b_max_lambda = 900
            else:
                raise RuntimeError('Invalid start frequency')

            freqs = freq0 + (np.arange(200) * 100e3 + 50e3)

            fit_data = np.loadtxt(fit_file)
            area_arcmin2 = fit_data[:, 4]
            cube_file = join('results', 'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
                             (freq0 / 1e6, freq1 / 1e6, weight),
                             'noise_l_cut_%i_%i_%s.fits'
                             % (b_min_lambda, b_max_lambda, weight))
            hdu_list = fits.open(cube_file)
            cube = hdu_list[0].data
            filename, ext = os.path.splitext(cube_file)
            hdu_list[0].header['BUNIT'] = 'K'
            hdu_list[0].header.set('comment', 'Converted to K from Jy/beam.')
            for i, freq in enumerate(freqs):
                image = cube[i, :, :]
                area_sr = (area_arcmin2[i] / 3600) * (pi / 180)**2
                image /= area_sr
                image *= (1e-26 * const.c.value**2) / (2 * const.k_B.value * freq**2)
            print('Writing cube:', filename + '_K.fits')
            hdu_list.writeto(filename + '_K.fits', clobber=True)


def main2():
    # http://goo.gl/Yh2RwR
    # http://goo.gl/ra2eWq

    beam_sigma = 50 * units.arcsec
    beam_area = 2 * np.pi * beam_sigma**2
    freq = 100 * units.MHz
    equiv = units.brightness_temperature(beam_area, freq)
    units.Jy.to(units.K, equivalencies=equiv)
    print('astropy answer: 1jy -> K = ', (1 * units.Jy).to(units.K, equivalencies=equiv))

    print('*' * 20)
    conv = const.c.value**2 / (2 * const.k_B.value * freq.to(units.Hz).value**2)
    print('beam area', beam_area, beam_area.to(units.sr))
    print('conv', conv)

    # 1 Jy = 1e-26 W (m^-2) (Hz^-2) -> K
    print('manual answer:', ((1.0 * 1.e-26 / beam_area.to(units.sr).value) * conv))

    print('*' * 20)
    beam_sigma = 40 * units.arcmin
    beam_area = 2 * np.pi * beam_sigma**2
    print(beam_area, beam_area.to(units.sr))
    print(beam_area.to(units.degree**2))
    # Convert area from arcmin2 to sr
    print((beam_area.value/3600) * (pi/180)**2)

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


if __name__ == '__main__':
    # main2()
    # print('--- ' * 5)
    main()

