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
from noise_cube.sim_cube import (element_effective_area, system_temp,
                                 evaluate_noise_rms,
                                 evaluate_scaled_noise_rms,
                                 raw_t_sys, raw_a_eff)
# seaborn.set_style('white')


def main():

    t_acc_scaled = (1000 * 3600) / 300
    t_acc = 5.0
    n = 100
    freqs = np.linspace(50e6, 650e6, n)
    t_sys, a_eff, sigma_pq = np.zeros(n), np.zeros(n), np.zeros(n)
    for i, freq in enumerate(freqs):
        t_sys[i] = system_temp(freq)
        a_eff[i] = element_effective_area(freq)
        sigma_pq[i] = evaluate_noise_rms(freq, t_acc=t_acc,
                                         bw_hz=100e3, eta=1.0,
                                         num_antennas=256)

    fig, (ax1, ax2) = plt.subplots(figsize=(10, 8), nrows=2, ncols=1,
                                        sharex=True)
    fig.subplots_adjust(left=0.1, right=0.96, bottom=0.08, top=0.96,
                        hspace=0.08, wspace=0)

    raw_freqs = np.array(raw_t_sys['freqs']) / 1e6
    ax1.plot(raw_freqs, raw_t_sys['values'], 'k+', ms=10, mew=2.0, label='data')
    ax1.plot(freqs / 1e6, t_sys, label='fit')
    t_sys_ = 60.0 * (const.c.value/freqs)**2.55
    ax1.plot(freqs / 1e6, t_sys_, 'r--', lw=2.0, label='60.0 $\lambda^{2.55}$')
    ax1.set_ylabel('System temperature (K)')
    ax1.legend(frameon=True)
    ax1.set_yscale('log')
    ax1.set_ylim(10)
    ax1.grid()

    ax2.plot(raw_freqs, raw_a_eff['values'], 'k+', ms=10, mew=2.0, label='data')
    ax2.plot(freqs / 1e6, a_eff, label='fit')
    a_eff_ = (const.c.value / freqs)**2 / 3
    ax2.plot(freqs / 1e6, a_eff_, 'r--', lw=2.0, label='$\lambda^2 / 3$')
    a_eff_ = np.ones_like(freqs) * math.pi * (1.5467851914 / 2)**2
    ax2.plot(freqs / 1e6, a_eff_, 'r:', lw=2.0, label='physical area (d=1.55m)')
    ax2.set_ylabel('Element effective area (m$^2$)')
    ax2.set_xlabel('Frequency (MHz)')
    ax2.legend(frameon=True)
    ax2.set_ylim(0.1, 2.2)
    ax2.grid()

    # ax3.plot(freqs / 1e6, sigma_pq)
    # ax3.set_ylabel('Baseline noise RMS $\sigma_{S}^{p,q}$ (Jy)')
    # ax3.grid()

    ax2.set_xlim(45, 355)

    fig.savefig('Tsys_Aeff.eps')

    plt.show()



    # element_effective_area(freq_hz)
    #
    # fig = plt.figure(figsize=(8, 8))
    # ax = fig.add_subplot(111)
    # ax.plot(freqs, a_eff, 'r.', ms=10.0, label='data')
    # freqs1_ = np.linspace(np.log10(0.05 * 1e9), np.log10(0.11 * 1e9), 500)
    # freqs2_ = np.linspace(np.log10(0.11 * 1e9), np.log10(0.65 * 1e9), 500)
    # ax.plot(10**freqs1_, 10**f1(freqs1_), 'b-', lw=1.0, label='interp')
    # ax.plot(10 ** freqs2_, 10 ** f2(freqs2_), 'b-', lw=1.0)
    # freqs_ = np.linspace(0.05 * 1e9, 0.65 * 1e9, 200)
    # areas = (const.c.value / freqs_)**2 / 3
    # ax.plot(freqs_, areas, 'g--', lw=1.0, label=r'$\lambda^{2} / 3$')
    # areas = np.ones_like(freqs_) * math.pi * (1.5 / 2) ** 2
    # ax.plot(freqs_, areas, 'g:', lw=1.5, label='area limit (d=1.5 m)')
    # ax.legend(fontsize='medium')
    # ax.set_xlabel('Frequency (hz)')
    # ax.set_ylabel(r'Element effective area ($m^2$)')
    # # ax.set_xscale('log')
    # # ax.set_yscale('log')
    # ax.set_ylim(0, 2)
    # ax.grid(True)
    # ax.set_xlim(freqs[0], freqs[-1])
    # plt.savefig('effective_area_lin.png')
    # plt.show()



    # element_effective_area(50e6, True)
    # system_temp(50e6, True)
    return

    freqs = 50e6 + np.arange(601) * 1e6
    s_pq_jy_1 = np.zeros(freqs.shape[0])
    s_pq_jy_2 = np.zeros(freqs.shape[0])
    for i, freq in enumerate(freqs):
        s_pq_jy_1[i] = evaluate_scaled_noise_rms(freq, 300, 100e3, 1.0, 1000)
        s_pq_jy_2[i] = evaluate_scaled_noise_rms(freq, 20, 100e3, 1.0, 1000)
    fig, ax = plt.subplots(figsize=(8, 8), nrows=1, ncols=1)
    ax.plot(freqs / 1e6, s_pq_jy_1, label='300 times')
    # Note: the 20 time version will be lower as this has a larger effective
    # t_acc
    # The resulting image noise from these two should be the same as the
    # 20 times will have 300/20 fewer points
    ax.plot(freqs / 1e6, s_pq_jy_2, label='20 times')
    ax.plot(freqs / 1e6, s_pq_jy_2 * 15**0.5, '--', label='20 times scaled')
    ax.set_xlim(0, 350)
    ax.set_xlabel('Frequency (MHz)')
    # TODO(BM) convert to expected image noise in K with simple beam area based
    #  on lambda / B
    ax.legend()
    ax.grid()
    plt.show()


if __name__ == '__main__':
    main()
