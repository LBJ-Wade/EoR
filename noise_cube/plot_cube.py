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
from scipy.interpolate import interp1d
import sys
from noise_cube.sim_cube import (evaluate_noise_rms,
                                 evaluate_scaled_noise_rms)
seaborn.set_style('ticks')
# seaborn.reset_orig()


def plot_noise(start_freq, l1, l2):
    fig, (ax1, ax2) = plt.subplots(figsize=(8, 8), ncols=1, nrows=2,
                                   sharex=True)

    # ---------------------------------
    weights = 'natural'
    algorithm = 'fft'
    results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    print(results_dir)
    noise_cube_file_k = join(results_dir, 'noise_%05.1fMHz_lcut_%04i_%04i_%s_%s_K.fits' % (start_freq/1e6, l1, l2, weights, algorithm))
    noise_cube_k = fits.getdata(noise_cube_file_k)
    noise_cube_std_k = np.std(noise_cube_k, axis=(1, 2))

    noise_cube_file = join(results_dir, 'noise_%05.1fMHz_lcut_%04i_%04i_%s_%s.fits' % (start_freq/1e6, l1, l2, weights, algorithm))
    noise_cube = fits.getdata(noise_cube_file)
    noise_cube_std = np.std(noise_cube, axis=(1, 2))

    noise_sigma_file = join(results_dir,
                            'noise_sigma_%05.1fMHz_lcut_%04i_%04i_%s_%s.npz' %
                            (start_freq/1e6, l1, l2, weights, algorithm))
    noise_sigma = np.load(noise_sigma_file)

    uvw_coords = np.load(join(results_dir, 'uvw_m.npz'))
    uu = uvw_coords['uu']
    num_times = uvw_coords['num_times']
    x = uvw_coords['x']
    print(x.shape[0])
    print(uu.shape[0])

    freqs = noise_sigma['freqs']
    sigma_im = noise_sigma['sigma_im']
    sigma_pq = noise_sigma['sigma_pq']
    num_coords = noise_sigma['num_coords']
    print('sigma_pq (Stokes-I) [Jy]', sigma_pq[:5])
    print('sigma_im (Stokes-I) [uJy/beam]',
          (sigma_pq[:5] / num_coords[:5]**0.5) * 1e6)

    # ax1.plot(freqs / 1e6, noise_cube_std_k * 1e3)
    ax1.plot(freqs / 1e6, noise_cube_std * 1e6, 'b-',
             label='$\sigma_{im}$ noise cube: %s, %s' % (algorithm, weights))
    ax1.plot(freqs / 1e6, sigma_im * 1e6, 'k--', label='$\sigma_{p,q}/\sqrt{n}$')

    # ---------------------------------
    # weights = 'natural'
    # algorithm = 'w-projection_64'
    # results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    # print(results_dir)
    # noise_cube_file_k = join(results_dir, 'noise_050.0MHz_lcut_0007_0251_%s_%s_K.fits' % (weights, algorithm))
    # noise_cube_k = fits.getdata(noise_cube_file_k)
    # noise_cube_std_k = np.std(noise_cube_k, axis=(1, 2))
    #
    # noise_cube_file = join(results_dir, 'noise_050.0MHz_lcut_0007_0251_%s_%s.fits' % (weights, algorithm))
    # noise_cube = fits.getdata(noise_cube_file)
    # noise_cube_std = np.std(noise_cube, axis=(1, 2))
    #
    # noise_sigma_file = join(results_dir,
    #                         'noise_sigma_%05.1fMHz_lcut_0007_0251_%s_%s.npz' %
    #                         (start_freq[0]/1e6, weights, algorithm))
    # noise_sigma = np.load(noise_sigma_file)
    #
    # uvw_coords = np.load(join(results_dir, 'uvw_m.npz'))
    # uu = uvw_coords['uu']
    # num_times = uvw_coords['num_times']
    # x = uvw_coords['x']
    # print(x.shape[0])
    # print(uu.shape[0])
    #
    # freqs = noise_sigma['freqs']
    # sigma_im = noise_sigma['sigma_im']
    # sigma_pq = noise_sigma['sigma_pq']
    # num_coords = noise_sigma['num_coords']
    # print('sigma_pq (Stokes-I) [Jy]', sigma_pq[:5])
    # print('sigma_im (Stokes-I) [uJy/beam]',
    #       (sigma_pq[:5] / num_coords[:5]**0.5) * 1e6)
    #
    # # ax1.plot(freqs / 1e6, noise_cube_std_k * 1e3)
    # ax1.plot(freqs / 1e6, noise_cube_std * 1e6, c='k', label='noise cube, %s' % algorithm)
    # ax1.plot(freqs / 1e6, sigma_im * 1e6, label='$\sigma_{im}$')
    # # ax1.plot(freqs / 1e6, (sigma_pq / num_coords**0.5) * 1e6, '--', c='g',
    # #          label=r'$\frac{\sigma_{p,q}}{\sqrt{n}}$')

    # ------------------
    weights = 'uniform'
    algorithm = 'fft'
    results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    print(results_dir)
    noise_cube_file_k = join(results_dir,
                             'noise_%05.1fMHz_lcut_%04i_%04i_%s_%s_K.fits' %
                             (start_freq/1e6, l1, l2, weights, algorithm))
    noise_cube_k = fits.getdata(noise_cube_file_k)
    noise_cube_std_k = np.std(noise_cube_k, axis=(1, 2))

    noise_cube_file = join(results_dir,
                           'noise_%05.1fMHz_lcut_%04i_%04i_%s_%s.fits' %
                           (start_freq / 1e6, l1, l2, weights, algorithm))
    noise_cube = fits.getdata(noise_cube_file)
    noise_cube_std = np.std(noise_cube, axis=(1, 2))

    noise_sigma_file = join(results_dir,
                            'noise_sigma_%05.1fMHz_lcut_%04i_%04i_%s_%s.npz' %
                            (start_freq / 1e6, l1, l2, weights, algorithm))
    noise_sigma = np.load(noise_sigma_file)

    uvw_coords = np.load(join(results_dir, 'uvw_m.npz'))
    uu = uvw_coords['uu']
    num_times = uvw_coords['num_times']
    x = uvw_coords['x']
    print(x.shape[0])
    print(uu.shape[0])

    freqs = noise_sigma['freqs']
    sigma_im = noise_sigma['sigma_im']
    sigma_pq = noise_sigma['sigma_pq']
    num_coords = noise_sigma['num_coords']
    print('sigma_pq (Stokes-I) [Jy]', sigma_pq[:5])
    print('sigma_im (Stokes-I) [uJy/beam]', (sigma_pq[:5] / num_coords[:5]**0.5) * 1e6)

    ax2.plot(freqs / 1e6, noise_cube_std * 1e6, 'r-',
             label='$\sigma_{im}$ noise cube: %s, %s' % (algorithm, weights))
    ax2.plot(freqs / 1e6, sigma_im * 1e6, 'k--', label='$\sigma_{p,q}/\sqrt{n}$')

    # # ------------------
    # weights = 'uniform'
    # algorithm = 'w-projection_64'
    # results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    # print(results_dir)
    # noise_cube_file_k = join(results_dir,
    #                          'noise_050.0MHz_lcut_0007_0251_%s_%s_K.fits' % (
    #                          weights, algorithm))
    # noise_cube_k = fits.getdata(noise_cube_file_k)
    # noise_cube_std_k = np.std(noise_cube_k, axis=(1, 2))
    #
    # noise_cube_file = join(results_dir,
    #                        'noise_050.0MHz_lcut_0007_0251_%s_%s.fits' % (
    #                        weights, algorithm))
    # noise_cube = fits.getdata(noise_cube_file)
    # noise_cube_std = np.std(noise_cube, axis=(1, 2))
    #
    # noise_sigma_file = join(results_dir,
    #                         'noise_sigma_%05.1fMHz_lcut_0007_0251_%s_%s.npz' %
    #                         (start_freq[0] / 1e6, weights, algorithm))
    # noise_sigma = np.load(noise_sigma_file)
    #
    # uvw_coords = np.load(join(results_dir, 'uvw_m.npz'))
    # uu = uvw_coords['uu']
    # num_times = uvw_coords['num_times']
    # x = uvw_coords['x']
    # print(x.shape[0])
    # print(uu.shape[0])
    #
    # freqs = noise_sigma['freqs']
    # sigma_im = noise_sigma['sigma_im']
    # sigma_pq = noise_sigma['sigma_pq']
    # num_coords = noise_sigma['num_coords']
    # print('sigma_pq (Stokes-I) [Jy]', sigma_pq[:5])
    # print('sigma_im (Stokes-I) [uJy/beam]', (sigma_pq[:5] / num_coords[:5]**0.5) * 1e6)
    #
    # ax2.plot(freqs / 1e6, noise_cube_std * 1e6, c='k',
    #          label='noise cube, %s' % algorithm)

    # # ax2.plot(freqs / 1e6, (sigma_pq / num_coords**0.5) * 1e6, '--',
    # #          label=r'$\frac{\sigma_{p,q}}{\sqrt{n}}$')


    # ---------------------------------
    # ax1.set_yscale('log')
    ax2_lim = ax2.get_ylim()
    print(ax2_lim)
    ax2.set_ylim(ax2_lim[0], ax2_lim[1]*1.2)
    ax1.grid()
    ax2.grid()
    ax1.legend(frameon=True)
    ax2.legend(frameon=True)
    ax1.set_title('Natural weights')
    ax2.set_title('Uniform weights')
    ax1.set_ylabel('Image RMS ($\mu$Jy/beam)')
    ax2.set_ylabel('Image RMS ($\mu$Jy/beam)')
    ax2.set_xlabel('Frequency (MHz)')
    fig.savefig('noise_%05.1fMHz.png' % (start_freq / 1e6))
    plt.show()


def plot_fit(start_freq):

    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(8, 10), nrows=3, sharex=True)
    ax1_y2 = ax1.twinx()
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.06, top=0.98,
                        hspace=0.1, wspace=0)

    # -----------------
    weights = 'uniform'
    algorithm = 'fft'
    results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    psf_fit_file = join(results_dir, 'psf_fit_%05.1fMHz_%s.npz' %
                            (start_freq / 1e6, weights))
    psf_fit = np.load(psf_fit_file)
    freqs = start_freq + np.arange(80) * 100e3 + 100e3 / 2
    theta_deg = psf_fit['theta_deg']
    sigma_x_arcmin = psf_fit['sigma_x_arcmin']
    sigma_y_arcmin = psf_fit['sigma_y_arcmin']
    area_arcmin2 = psf_fit['area_arcmin2']
    fit_rms = psf_fit['fit_rms']

    l1 = ax1.plot(freqs / 1e6, sigma_x_arcmin, '-', c='r', label='%s weights, $\sigma_{x}$' % weights)
    l2 = ax1.plot(freqs / 1e6, sigma_y_arcmin, '--', c='r', label='%s weights, $\sigma_{y}$' % weights)
    l3 = ax1_y2.plot(freqs / 1e6, theta_deg, ':', color='r', label=r'%s weights, $\theta$' % weights)
    ax2.plot(freqs / 1e6, area_arcmin2, c='r', label='%s weights' % weights)
    # ax3.plot(freqs / 1e6, fit_rms, c='r', label='%s %s' % (weights, algorithm))

    # -----------------
    weights = 'natural'
    algorithm = 'fft'
    results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    psf_fit_file = join(results_dir, 'psf_fit_%05.1fMHz_%s.npz' %
                            (start_freq / 1e6, weights))
    psf_fit = np.load(psf_fit_file)
    freqs = start_freq + np.arange(80) * 100e3 + 100e3 / 2
    theta_deg = psf_fit['theta_deg']
    sigma_x_arcmin = psf_fit['sigma_x_arcmin']
    sigma_y_arcmin = psf_fit['sigma_y_arcmin']
    area_arcmin2 = psf_fit['area_arcmin2']
    fit_rms = psf_fit['fit_rms']

    l4 = ax1.plot(freqs / 1e6, sigma_x_arcmin, '-', c='b', label='%s weights, $\sigma_{x}$' % weights)
    l5 = ax1.plot(freqs / 1e6, sigma_y_arcmin, '--', c='b', label='%s weights, $\sigma_{y}$' % weights)
    l6 = ax1_y2.plot(freqs / 1e6, theta_deg, ':', color='b', label=r'%s weights, $\theta$' % weights)
    ax3.plot(freqs / 1e6, area_arcmin2, c='b', label='%s weights' % weights)
    # ax3.plot(freqs / 1e6, fit_rms, c='b', label='%s %s' % (weights, algorithm))


    # # -----------------
    # weights = 'uniform'
    # algorithm = 'fft'
    # results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    # psf_fit_file = join(results_dir, 'psf_fit_%05.1fMHz_%s.npz' %
    #                         (start_freq / 1e6, weights))
    # psf_fit = np.load(psf_fit_file)
    # freqs = start_freq[0] + np.arange(80) * 100e3 + 100e3 / 2
    # theta_deg = psf_fit['theta_deg']
    # sigma_x_arcmin = psf_fit['sigma_x_arcmin']
    # sigma_y_arcmin = psf_fit['sigma_y_arcmin']
    # area_arcmin2 = psf_fit['area_arcmin2']
    # fit_rms = psf_fit['fit_rms']
    #
    # ax2.plot(freqs / 1e6, area_arcmin2, '--', c='k', label='%s %s' % (weights, algorithm))
    # ax3.plot(freqs / 1e6, fit_rms, '--', c='k', label='%s %s' % (weights, algorithm))
    #
    # # -----------------
    # weights = 'natural'
    # algorithm = 'fft'
    # results_dir = 'noise_cubes_%s_%s' % (weights, algorithm)
    # psf_fit_file = join(results_dir, 'psf_fit_%05.1fMHz_%s.npz' %
    #                         (start_freq[0] / 1e6, weights))
    # psf_fit = np.load(psf_fit_file)
    # freqs = start_freq[0] + np.arange(80) * 100e3 + 100e3 / 2
    # theta_deg = psf_fit['theta_deg']
    # sigma_x_arcmin = psf_fit['sigma_x_arcmin']
    # sigma_y_arcmin = psf_fit['sigma_y_arcmin']
    # area_arcmin2 = psf_fit['area_arcmin2']
    # fit_rms = psf_fit['fit_rms']
    #
    # ax2.plot(freqs / 1e6, area_arcmin2, '--', c='r', label='%s %s' % (weights, algorithm))
    # ax3.plot(freqs / 1e6, fit_rms, '--', c='r', label='%s %s' % (weights, algorithm))
    #

    # -----------------
    #ax1.grid()
    ax1_lim = ax1.get_ylim()
    ax1.set_ylim(ax1_lim[0], ax1_lim[1] + math.fabs(ax1_lim[1]*0.4))
    ax1_y2_lim = ax1_y2.get_ylim()
    ax1_y2.set_ylim(ax1_y2_lim[0], ax1_y2_lim[1] + math.fabs(ax1_y2_lim[1]*0.2))
    ax1.xaxis.grid()
    ax2.grid()
    ax3.grid()
    lns = l1 + l4 + l2 + l5 + l3 + l6
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, bbox_to_anchor=(0.5, 1.0),
               loc='upper center',
               ncol=3,
               frameon=True,
               fontsize='small')

    # ax1.legend(lns, labs, frameon=True)
    ax2.legend(frameon=True, fontsize='small')
    ax3.legend(frameon=True, fontsize='small')
    # ax3.legend(frameon=True, loc='center right', fontsize='small')
    ax1_y2.set_ylabel('PSF main lobe fit angle (degrees)')
    ax1.set_ylabel('PSF main lobe fit $\sigma$ (arcmin)')
    ax2.set_ylabel('PSF main lobe area (arcmin$^{2}$)')
    ax3.set_ylabel('PSF main lobe area (arcmin$^{2}$)')
    # ax3.set_ylabel('Fit RMS')
    ax3.set_xlabel('Frequency (MHz)')
    fig.savefig('beam_fit_%05.1fMHz.png' % (start_freq / 1e6))
    plt.show()


if __name__ == '__main__':
    # freq0 = 50e6
    # l1, l2 = 7, 251
    # freq0 = 98e6
    # l1, l2 = 13, 493
    # freq0 = 146e6
    # l1, l2 = 18, 734
    freq0 = 194e6
    l1, l2 = 24, 975

    plot_fit(start_freq=freq0)
    plot_noise(start_freq=freq0, l1=l1, l2=l2)
