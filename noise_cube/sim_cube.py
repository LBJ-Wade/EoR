# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
import oskar
from os.path import join
import os
from astropy import constants as const
import numpy as np
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn
from pyuvwsim import convert_enu_to_ecef, evaluate_baseline_uvw
import time
import sys
from astropy.io import fits
seaborn.set_style('ticks')


def element_effective_area(freq_hz, debug_plot=False):
    """Return SKA1 Low element pattern effective area for given frequency

    Effective area values provided by Eloy de Lera Acedo
    (email: eloy .at. mrao.cam.ac.uk)

    Args:
        freq_hz (float): Frequency, in Hz

    Returns:
        Element effective area, in m^2
    """
    freqs = [0.05, 0.07, 0.11, 0.17, 0.25, 0.35, 0.45, 0.55, 0.65]
    a_eff = [1.8791, 1.8791, 1.8694, 1.3193, 0.6080, 0.2956, 0.2046,
             0.1384, 0.0792]
    freqs, a_eff = np.array(freqs) * 1e9, np.array(a_eff)
    f_cut = 2
    f1 = interp1d(np.log10(freqs[:f_cut+1]), np.log10(a_eff[:f_cut+1]),
                  kind='slinear')
    f2 = interp1d(np.log10(freqs[f_cut:]), np.log10(a_eff[f_cut:]),
                  kind='cubic')
    if debug_plot:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        ax.plot(freqs, a_eff, 'r.', ms=10.0, label='data')
        freqs1_ = np.linspace(np.log10(0.05 * 1e9), np.log10(0.11 * 1e9), 500)
        freqs2_ = np.linspace(np.log10(0.11 * 1e9), np.log10(0.65 * 1e9), 500)
        ax.plot(10**freqs1_, 10**f1(freqs1_), 'b-', lw=1.0, label='interp')
        ax.plot(10 ** freqs2_, 10 ** f2(freqs2_), 'b-', lw=1.0)
        freqs_ = np.linspace(0.05 * 1e9, 0.65 * 1e9, 200)
        areas = (const.c.value / freqs_)**2 / 3
        ax.plot(freqs_, areas, 'g--', lw=1.0, label=r'$\lambda^{2} / 3$')
        areas = np.ones_like(freqs_) * math.pi * (1.5 / 2) ** 2
        ax.plot(freqs_, areas, 'g:', lw=1.5, label='area limit (d=1.5 m)')
        ax.legend(fontsize='medium')
        ax.set_xlabel('Frequency (hz)')
        ax.set_ylabel(r'Element effective area ($m^2$)')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0, 2)
        ax.grid(True)
        ax.set_xlim(freqs[0], freqs[-1])
        plt.show()
    if freq_hz <= freqs[f_cut]:
        return 10**f1(np.log10(freq_hz))
    else:
        return 10**f2(np.log10(freq_hz))


def system_temp(freq_hz, debug_plot=False):
    """Return SKA1 Low system temperatures for a given frequency.

    Values provided by Eloy de Lera Acedo
    (email: eloy .at. mrao.cam.ac.uk)

    Args:
        freq_hz (float): Frequency, in Hz

    Returns:
        System temperature, in K
    """
    freqs = [0.05, 0.07, 0.11, 0.17, 0.25, 0.35, 0.45, 0.55, 0.65]
    t_sys = [4.0409e3, 1.5029e3, 0.6676e3, 0.2936e3, 0.1402e3, 0.0873e3,
             0.0689e3, 0.0607e3, 0.0613e3]
    freqs, t_sys = np.array(freqs) * 1e9, np.array(t_sys)
    # f = interp1d(freqs, t_sys, kind='slinear')
    f2 = interp1d(np.log10(freqs), np.log10(t_sys), kind='cubic')
    if debug_plot:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        ax.plot(freqs, t_sys, 'r.', ms=10.0, label='data')
        #freqs_ = np.linspace(0.05 * 1e9, 0.65 * 1e9, 500)
        #ax.plot(freqs_, f(freqs_), 'b-', lw=1.0, label='interp')
        freqs_ = np.linspace(np.log10(0.05 * 1e9), np.log10(0.65 * 1e9), 500)
        ax.plot(10**freqs_, 10**f2(freqs_), 'b-', lw=1.0, label='interp')
        freqs_ = np.linspace(0.05 * 1e9, 0.65 * 1e9, 20)
        t_sys_ = 60.0 * (const.c / freqs_)**2.55
        ax.plot(freqs_, t_sys_, 'g--', lw=1.0, ms=5.0, mew=1.2,
                label='$60\lambda^{2.55}$')

        ax.legend(fontsize='large')
        ax.set_xlabel('Frequency (hz)')
        ax.set_ylabel('System temperature (K)')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_xlim(freqs[0], freqs[-1])
        plt.show()
    return 10**f2(np.log10(freq_hz))


def plot_sigma_im():
    eta = 1.0
    bw_hz = 100e3
    num_antennas = 256
    t_acc = 5.0
    obs_length_h = 1000.0
    num_times = (obs_length_h * 3600.0) / t_acc
    num_stations = 200
    num_baselines = (num_stations * (num_stations - 1)) // 2
    n = 2 * num_times * num_baselines
    freqs = np.linspace(50.0e6, 350.0e6, 2**16 // 128)
    sigma_pq = np.zeros_like(freqs)
    sigma_im = np.zeros_like(freqs)
    for i, freq in enumerate(freqs):
        t_sys = system_temp(freq)
        a_eff = element_effective_area(freq) * num_antennas
        sefd = (2.0 * const.k_B.value * t_sys * eta) / a_eff
        sefd *= 1e26  # Convert to Jy
        sigma_pq[i] = sefd / (2.0 * bw_hz * t_acc) ** 0.5
        sigma_im[i] = (sigma_pq[i] / (2 * n)**0.5) * 1e6
    print(sigma_im.min(), sigma_im.max())
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(freqs, sigma_im)
    ax.set_ylabel('Expected image RMS ($\mu$Jy/beam)')
    ax.set_xlabel('Frequency (Hz)')
    ax.grid(True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(sigma_im.min()*0.9, sigma_im.max()*1.1)
    ax.set_xlim(freqs[0], freqs[-1])
    plt.show()


def evaluate_noise_rms(freq_hz, t_acc=5.0, bw_hz=100e3, eta=1.0,
                       num_stations=200, obs_length_h=1000.0, verbose=False):
    debug_plot = False
    station_d = 35.0  # m
    num_antennas = 256  # Number of antennas per station

    t_sys = system_temp(freq_hz, debug_plot)
    a_eff = element_effective_area(freq_hz, debug_plot) * num_antennas

    # Single receiver polarisation SEFD
    sefd = (2.0 * const.k_B.value * t_sys * eta) / a_eff
    sefd *= 1e26  # Convert to Jy
    sigma_pq = sefd / (2.0 * bw_hz * t_acc) ** 0.5

    # Expected image noise (Stokes-I)
    num_times = (obs_length_h * 3600.0) / t_acc
    num_baselines = (num_stations * (num_stations - 1)) // 2
    n = 2 * num_times * num_baselines
    sigma_im = sigma_pq / n**0.5

    if verbose:
        print('*' * 60)
        print('- freq = %.2f MHz' % (freq_hz / 1e6))
        print('- t_sys = %.2f K' % t_sys)
        print('- station a_eff = %.2f m^2' % a_eff)
        print('- station area = %.2f m^2' % (math.pi * (station_d / 2) ** 2))
        print('- bandwidth = %.2f kHz' % (bw_hz / 1e3))
        print('- dump time = %.2f s' % t_acc)
        print('- SEFD = %.3f Jy' % sefd)
        print('- sigma_pq = %.3f Jy' % sigma_pq)
        print('- no. times (%.1f h) = %i' % (obs_length_h, int(num_times)))
        print('- no. stations = %i' % (int(num_stations)))
        print('- no. baselines = %i' % (int(num_baselines)))
        print('- expected %.1fh image noise = %.3f uJy' %
              (obs_length_h, sigma_im * 1e6))
        print('*' * 60)

    return sigma_pq, sigma_im


def evaluate_effective_noise_rms(
        freq_hz, num_times=48, bw_hz=100e3, eta=1.0, num_stations=200,
        obs_length_h=6.0, target_obs_length_h=1000.0, verbose=False):

    debug_plot = False
    station_d = 35.0  # m
    num_antennas = 256  # Number of antennas per station
    t_acc = (target_obs_length_h * 3600) / num_times
    t_sys = system_temp(freq_hz, debug_plot)
    a_eff = element_effective_area(freq_hz, debug_plot) * num_antennas

    # Single receiver polarisation SEFD
    sefd = (2.0 * const.k_B.value * t_sys * eta) / a_eff
    sefd *= 1e26  # Convert to Jy
    sigma_pq = sefd / (2.0 * bw_hz * t_acc) ** 0.5

    # Expected image noise (Stokes-I)
    num_baselines = (num_stations * (num_stations - 1)) // 2
    n = 2 * num_times * num_baselines
    sigma_im = sigma_pq / n**0.5

    if verbose:
        print('*' * 60)
        print('- freq = %.2f MHz' % (freq_hz / 1e6))
        print('- t_sys = %.2f K' % t_sys)
        print('- station a_eff = %.2f m^2' % a_eff)
        print('- station area = %.2f m^2' % (math.pi * (station_d / 2) ** 2))
        print('- bandwidth = %.2f kHz' % (bw_hz / 1e3))
        print('- dump time = %.2f s' % t_acc)
        print('- SEFD = %.3f Jy' % sefd)
        print('- sigma_pq = %.3f Jy' % sigma_pq)
        print('- no. times (%.1f h) = %i' % (obs_length_h, int(num_times)))
        print('- no. stations = %i' % (int(num_stations)))
        print('- no. baselines = %i' % (int(num_baselines)))
        print('- expected %.1fh image noise = %.3f uJy' %
              (target_obs_length_h, sigma_im * 1e6))
        print('*' * 60)

    return sigma_pq, sigma_im


def plot_telescope(x, y, r_cut, station_d=35.0):
    # Plot the telescope model
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    c = plt.Circle((0, 0), r_cut, fill=False, color='r')
    ax.add_artist(c)
    for x_, y_ in zip(x, y):
        c = plt.Circle((x_, y_), station_d/2, fill=False)
        ax.add_artist(c)
    ax.grid()
    ax.set_xlim(-r_cut*1.2, r_cut*1.2)
    ax.set_ylim(-r_cut*1.2, r_cut*1.2)
    ax.set_xlabel('east (m)')
    ax.set_ylabel('north (m)')
    fig.savefig('stations.png')
    plt.close(fig)


def plot_uv_coords(uu, vv, cut_r=None, units='m', filename=None):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(uu, vv, 'k.', ms=2.0, alpha=0.01)
    ax.plot(-uu, -vv, 'k.', ms=2.0, alpha=0.01)
    if cut_r:
        c = plt.Circle((0, 0), cut_r, fill=False, color='r', lw=2.0, alpha=1.0)
        ax.add_artist(c)
    ax.grid()
    ax.set_xlabel('uu (%s)' % units)
    ax.set_ylabel('vv (%s)' % units)
    if not filename:
        filename = 'uv_scatter.png'
    fig.savefig(filename)
    plt.close(fig)


def load_telescope(r_cut, station_d, lon, lat, plot=False):
    # Load telescope model
    coords = np.loadtxt(join('v5.tm', 'layout.txt'))
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r = (x**2 + y**2)**0.5

    x = x[np.where(r < r_cut)]
    y = y[np.where(r < r_cut)]
    z = z[np.where(r < r_cut)]
    num_stations = x.shape[0]
    if plot:
        plot_telescope(x, y, r_cut, station_d=station_d)
    x, y, z = convert_enu_to_ecef(x, y, z, lon, lat)
    return x, y, z


def gen_uvw_m(x, y, z, obs_length_h, mjd_mid, t_acc, ra, dec, plot=False):
    num_times = int((obs_length_h * 3600.0) / t_acc)
    mjd_start = mjd_mid - ((obs_length_h / 2) / 24)
    uu, vv, ww = np.array([]), np.array([]), np.array([])
    print('- Simulating uv coordinates for %i times over %.1f h'
          % (num_times, obs_length_h))
    t0 = time.time()
    for i in range(num_times):
        mjd = mjd_start + (i * t_acc + t_acc / 2) / 86400
        uu_, vv_, ww_ = evaluate_baseline_uvw(x, y, z, ra, dec, mjd)
        uu = np.concatenate((uu, uu_))
        vv = np.concatenate((vv, vv_))
        ww = np.concatenate((ww, ww_))
    print('- Generated uv %i coordinates in %.2f s'
          % (uu.shape[0], time.time() - t0))
    if plot:
        plot_uv_coords(uu, vv)
    return uu, vv, ww, num_times


def main():
    # Simulate 1000h observation as a 4-8h observation with 1000h noise.
    # TODO(BM) compare t_sys with GSM value?
    # TODO(BM) Define observation track ra, dec, mjd / ha range.
    freqs_hz = np.linspace(100e6, 120e6, 10)
    # freqs_hz = np.array([50e6])
    lon, lat = 116.63128900, -26.69702400
    ra, dec = 68.698903779331502, -26.568851215532160
    mjd_mid = 57443.4375000000
    station_d = 35.0
    r_cut_uv_l = None  # cut radius in wavelengths or None to disable
    # TODO(BM) add an inner uv cut as well as outer?
    # TODO(BM) dont want to do a lambda cut across the whole band but
    #          for channels around certain freqs as in Cath T's bandpass memo
    obs_length_h = 6.0
    target_obs_length_h = 1000.0
    t_acc = 300.0
    bw_hz = 100e3
    im_size = 1024  # TODO(BM) set this properly!
    fov_deg = 10.0  # TODO(BM) set this properly!
    if not os.path.isdir('uv_plots'):
        os.makedirs('uv_plots')

    # hpbw = np.degrees(wavelength / station_d)
    # size = (r_cut * 2) / station_d

    # Load the telescope model.
    r_cut = r_cut_uv_l * (const.c.value / freqs_hz.min()) \
        if r_cut_uv_l else 1.5e3
    x, y, z = load_telescope(r_cut, station_d, lon, lat, plot=True)

    # Generate uvw coordinates in m
    uu, vv, ww, num_times = gen_uvw_m(x, y, z, obs_length_h, mjd_mid, t_acc,
                                      ra, dec, plot=True)

    # evaluate_noise_rms(freqs_hz[-1], num_stations=x.shape[0], verbose=True)
    # evaluate_effective_noise_rms(freqs_hz[-1], num_stations=x.shape[0],
    #                              obs_length_h=obs_length_h,
    #                              num_times=num_times,
    #                              verbose=True)

    # Create imager object
    imager = oskar.imager.Imager('single')
    images = np.empty((freqs_hz.shape[0], im_size, im_size))

    for i, freq_hz in enumerate(freqs_hz):
        t0 = time.time()
        print('%3i - %8.2f MHz' % (i, freq_hz / 1e6), end=' ')

        # Scale to wavelengths and lambda cut
        wavelength = const.c.value / freq_hz
        uu_l = np.copy(uu) / wavelength
        vv_l = np.copy(vv) / wavelength
        ww_l = np.copy(ww) / wavelength
        if r_cut_uv_l:
            r_uv = (uu_l**2 + vv_l**2)**0.5
            uu_l = uu_l[np.where(r_uv < r_cut_uv_l)]
            vv_l = vv_l[np.where(r_uv < r_cut_uv_l)]
            ww_l = ww_l[np.where(r_uv < r_cut_uv_l)]
            plot_uv_coords(uu_l * wavelength, vv_l * wavelength,
                           r_cut_uv_l * wavelength,
                           filename=join('uv_plots',
                                         'uv_scatter_%06.2fMHz.png'
                                         % (freq_hz/1e6)))

        # Generate noise visibility amplitudes
        sigma_pq, sigma_im = evaluate_effective_noise_rms(
            freq_hz, num_times=num_times, bw_hz=bw_hz, num_stations=x.shape[0],
            obs_length_h=obs_length_h, target_obs_length_h=target_obs_length_h,
            verbose=False)
        # sigma_pq /= 2**0.5
        # amp = (np.random.randn(uu_l.shape[0]) * sigma_pq +
        #        1.0j * np.random.randn(uu_l.shape[0]) * sigma_pq)
        amp = np.ones(uu_l.shape[0], dtype='c16')

        image = imager.make_image(uu_l, vv_l, ww_l, amp, np.ones_like(uu_l),
                                  fov_deg, im_size)
        images[i, :, :] = image

        print(': image noise ~ %6.1f uJy/beam' % (sigma_im * 1e6), end=' ')
        print('(%6.1f)' % (np.std(image) * 1e6), end=' ')
        print(': %.2f s' % (time.time() - t0))
        sys.stdout.flush()

    # TODO convert units jy/beam -> K
    t0 = time.time()
    hdu = fits.PrimaryHDU(images)
    hdulist = fits.HDUList([hdu])
    if r_cut_uv_l:
        hdulist.writeto('test_l_cut_%.1f.fits' % r_cut_uv_l, clobber=True)
    else:
        hdulist.writeto('test_no_l_cut.fits', clobber=True)
    print('write fits: %.2f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
    # plot_sigma_im()
    # element_effective_area(50.0e6, debug_plot=True)
    # system_temp(50e6, debug_plot=True)



