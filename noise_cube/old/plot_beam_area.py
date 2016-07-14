# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from os.path import join, dirname
import seaborn
from math import floor, ceil, degrees, pi, log, sin, asin
from astropy import constants as const
from scipy.signal import savgol_filter
seaborn.set_style('ticks')


def smooth(x, window_len=11, window='hanning'):
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')
    y = np.convolve(w/w.sum(), s, mode='valid')
    return y


def main():
    window_length = 41
    # weight = 'uniform'
    weight = 'natural'
    freq0 = 200e6
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

    fit_data = np.loadtxt(fit_file)
    area = fit_data[:, 4]
    freq = freq0 + (np.arange(200) * 100e3 + 50e3)
    sigma_ = degrees(1 / b_max_lambda) * 60 / (2*(2*log(2))**0.5)
    area_ = 2 * pi * sigma_**2
    print(sigma_, fit_data[:5, 0])
    print(area_, area[:5])
    fig, ax = plt.subplots(figsize=(8, 8), ncols=1, nrows=1)
    ax.plot(freq / 1e6, area, '.', label='data')
    y = smooth(area, window_length, window='hanning')
    ax.plot(freq / 1e6,
            y[(y.shape[0] - freq.shape[0])//2:-(y.shape[0] - freq.shape[0])//2],
            '-', label='hanning smoothed')
    # https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
    y2 = savgol_filter(area, window_length, 1)
    ax.plot(freq / 1e6, y2, '-', label='Savgol filter')

    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Beam area (arcmin$^2$)')
    # ax.set_ylim(floor(area.min()), ceil(area.max()))
    ax.set_ylim(area.min(), area.max())
    ax.grid(True)
    ax.legend()
    fig.savefig(join(dirname(fit_file), 'beam_area_plot.png'))
    plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()

