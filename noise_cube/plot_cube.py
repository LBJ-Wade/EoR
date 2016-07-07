# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from os.path import join
import os
from astropy import constants as const
import numpy as np
import math
import matplotlib.pyplot as plt
import time
import sys
import shutil
from astropy.io import fits
from progressbar import ProgressBar, Percentage, Bar
import seaborn
seaborn.set_style('ticks')


if __name__ == '__main__':
    cube_u = fits.getdata(join('results',
                               'noise_50-70MHz_100kHz_5h_60s_uniform',
                               'noise_l_cut_10_230_uniform.fits'))
    cube_n = fits.getdata(join('results',
                               'noise_50-70MHz_100kHz_5h_60s_natural',
                               'noise_l_cut_10_230_natural.fits'))
    sigma_im_n = np.loadtxt(join('results',
                                 'noise_50-70MHz_100kHz_5h_60s_natural',
                                 'sigma_im.txt'))

    freq = 50e6 + (np.arange(200) * 100e3 + 50e3)
    std_u = np.zeros(cube_u.shape[0])
    std_n = np.zeros(cube_n.shape[0])
    pbar = ProgressBar(maxval=200, widgets=[Bar(), Percentage()])
    pbar.start()
    for i in range(cube_u.shape[0]):
        std_u[i] = np.std(cube_u[i, :, :])
        std_n[i] = np.std(cube_n[i, :, :])
        time.sleep(0.001)
        pbar.update(i)
    pbar.finish()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(freq / 1e6, std_u * 1e6, 'r-', label='uniform')
    ax.plot(freq / 1e6, std_n * 1e6, 'k-', label='natural')
    ax.plot(freq / 1e6, sigma_im_n * 1e6, 'g-',
            label=r'$\sigma_{\mathrm{image}}$')
    ax.set_ylabel(r'Image rms ($\mu$Jy/beam)')
    ax.set_xlabel('Frequency (MHz)')
    ax.legend(fontsize='medium')
    ax.grid()
    plt.show()
