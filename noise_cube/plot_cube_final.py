# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals
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


def plot_noise():
    start_freq = [50e6]

    results_dir = 'results_noise_cubes_uniform_w-projection'
    noise_cube_file_k = join(results_dir, 'noise_050.0MHz_lcut_0009_0251_uniform_w-projection_K.fits')
    noise_cube_k = fits.getdata(noise_cube_file_k)
    noise_cube_std_k = np.std(noise_cube_k, axis=(1, 2))

    noise_sigma_file = join(results_dir, 'noise_sigma_%05.1fMHz.npz' %
                            (start_freq[0]/1e6))
    noise_sigma = np.load(noise_sigma_file)
    freqs = noise_sigma['freqs']
    sigma_im = noise_sigma['sigma_im']

    fig, ax1 = plt.subplots(figsize=(8, 8))
    ax1.plot(freqs / 1e6, noise_cube_std_k * 1e3)
    # ax1.set_yscale('log')
    ax1.set_ylabel('Image rms (mK)')
    plt.show()


def plot_fit():
    start_freq = [50e6]
    results_dir = 'results_noise_cubes_uniform_w-projection'
    psf_fit_file = join(results_dir, 'psf_fit_%05.1fMHz_uniform.npz' %
                            (start_freq[0]/1e6))
    psf_fit = np.load(psf_fit_file)

    freqs = start_freq[0] + np.arange(200) * 100e3 + 100e3 / 2
    theta_deg = psf_fit['theta_deg']
    sigma_x_arcmin = psf_fit['sigma_x_arcmin']
    sigma_y_arcmin = psf_fit['sigma_y_arcmin']
    area_arcmin2 = psf_fit['area_arcmin2']
    fit_rms = psf_fit['fit_rms']

    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(12, 12), nrows=3)
    l1 = ax1.plot(freqs / 1e6, sigma_x_arcmin, label='$\sigma_{x}$')
    l2 = ax1.plot(freqs / 1e6, sigma_y_arcmin, label='$\sigma_{y}$')

    ax1_y2 = ax1.twinx()
    l3 = ax1_y2.plot(freqs / 1e6, theta_deg, ':', color='r', label='theta')

    lns = l1 + l2 + l3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0)

    ax1.set_ylabel('beam sigma (arcmin)')
    ax1_y2.set_ylabel('beam angle (degrees)')

    ax2.plot(freqs / 1e6, area_arcmin2)
    ax2.set_ylabel('beam area (arcmin$^{2}$)')

    ax3.plot(freqs / 1e6, fit_rms)
    ax3.set_ylabel('FIT RMS')

    plt.show()


if __name__ == '__main__':
    plot_fit()
    plot_noise()
