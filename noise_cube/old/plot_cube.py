# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from os.path import join
import os
from astropy import constants as const
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from astropy.io import fits
from progressbar import ProgressBar, Percentage, Bar
import seaborn
import sys
seaborn.set_style('ticks')


def get_fits_file(out_dir, f0, f1, weight):
    dir = join(out_dir, 'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
               (f0, f1, weight))
    file = [f for f in os.listdir(dir) if 'noise_' in f]
    return join(dir, file[0])


def get_sigma_im_file(out_dir, f0, f1, weight):
    dir = join(out_dir, 'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
               (f0, f1, weight))
    file = [f for f in os.listdir(dir) if 'sigma_im.txt' in f]
    return join(dir, file[0])


def main2():
    for freq0 in (50e6, 120e6, 200e6):
        freqs = freq0 + (np.arange(200) * 100e3 + 50e3)
        freq1 = freq0 + 20e6
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

        fig, ax = plt.subplots(figsize=(8, 8))
        # for weight in ('uniform', 'natural'):
        for weight in ['natural']:
            print(freq0/1e6,  weight)
            # cube_file_jy = join('results',
            #                     'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
            #                     (freq0 / 1e6, freq1 / 1e6, weight),
            #                     'noise_l_cut_%i_%i_%s.fits'
            #                     % (b_min_lambda, b_max_lambda, weight))
            # cube_jy = fits.getdata(cube_file_jy)
            cube_file_k = join('results',
                               'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
                               (freq0 / 1e6, freq1 / 1e6, weight),
                               'noise_l_cut_%i_%i_%s_K.fits'
                               % (b_min_lambda, b_max_lambda, weight))
            cube_k = fits.getdata(cube_file_k)
            img_std = np.std(cube_k, axis=(1, 2))
            ax.plot(freqs / 1e6, img_std * 1e3, label=weight)
        ax.set_ylabel(r'Image rms (mK)')
        ax.set_xlabel('Frequency (MHz)')
        ax.legend(fontsize='medium')
        ax.grid()
        fig.savefig(join('results', 'noise_%iMHz_K_n.png' % (freq0/1e6)))
        plt.close(fig)


def main():
    out_dir = 'results'
    f0, f1 = 200, 220
    cube_u = fits.getdata(get_fits_file(out_dir, f0, f1, 'uniform'))
    cube_n = fits.getdata(get_fits_file(out_dir, f0, f1, 'natural'))
    sigma_im = np.loadtxt(get_sigma_im_file(out_dir, f0, f1, 'natural'))

    #TODO(BM) get from fits header
    freq = f0 * 1e6 + (np.arange(200) * 100e3 + 50e3)
    std_u = np.zeros(cube_u.shape[0])
    std_n = np.zeros(cube_n.shape[0])
    pbar = ProgressBar(maxval=200, widgets=[Bar(), Percentage()])
    pbar.start()
    for i in range(cube_u.shape[0]):
        std_u[i] = np.std(cube_u[i, :, :])
        std_n[i] = np.std(cube_n[i, :, :])
        pbar.update(i)
    pbar.finish()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(freq / 1e6, std_u * 1e6, 'r-', label='uniform')
    ax.plot(freq / 1e6, std_n * 1e6, 'k-', label='natural')
    ax.plot(freq / 1e6, sigma_im * 1e6, 'g-',
            label=r'$\sigma_{\mathrm{image}}$')
    ax.set_ylabel(r'Image rms ($\mu$Jy/beam)')
    ax.set_xlabel('Frequency (MHz)')
    ax.legend(fontsize='medium')
    ax.grid()
    fig.savefig(join('results', 'noise_%iMHz.png' % f0))
    plt.close(fig)


if __name__ == '__main__':
    # main()
    main2()
