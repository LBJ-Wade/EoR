# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from os.path import join
import os
import math
import numpy as np
from numpy.random import randn
from math import ceil, floor, sin, cos, asin, degrees, radians, pi, log
from scipy.interpolate import interp1d
from scipy import optimize as opt
import time
from math import sin, degrees, radians, asin
from astropy import constants as const
from astropy.io import fits
import ephem
from pyuvwsim import convert_enu_to_ecef, evaluate_baseline_uvw
from oskar.imager import Imager
import seaborn
import matplotlib.pyplot as plt
from progressbar import ProgressBar, Percentage, Bar, ETA, Timer, Counter
seaborn.set_style('ticks')


raw_a_eff = {
    'freqs': [0.05e9, 0.07e9, 0.11e9, 0.17e9, 0.25e9, 0.35e9, 0.45e9,
              0.55e9, 0.65e9],
    'values': [1.8791, 1.8791, 1.8694, 1.3193, 0.6080, 0.2956, 0.2046,
               0.1384, 0.0792]
}
raw_t_sys = {
    'freqs': [0.05e9, 0.07e9, 0.11e9, 0.17e9, 0.25e9, 0.35e9, 0.45e9,
              0.55e9, 0.65e9],
    'values': [4.0409e3, 1.5029e3, 0.6676e3, 0.2936e3, 0.1402e3, 0.0873e3,
               0.0689e3, 0.0607e3, 0.0613e3]
}


def element_effective_area(freq_hz):
    """Return SKA1 Low element pattern effective area for given frequency

    Effective area values provided by Eloy de Lera Acedo
    (email: eloy .at. mrao.cam.ac.uk)

    Args:
        freq_hz (float): Frequency, in Hz

    Returns:
        Element effective area, in m^2
    """
    freqs, a_eff = np.array(raw_a_eff['freqs']), np.array(raw_a_eff['values'])
    f_cut = 2
    f1 = interp1d(np.log10(freqs[:f_cut+1]), np.log10(a_eff[:f_cut+1]),
                  kind='slinear')
    f2 = interp1d(np.log10(freqs[f_cut:]), np.log10(a_eff[f_cut:]),
                  kind='cubic')
    if freq_hz <= freqs[f_cut]:
        return 10**f1(np.log10(freq_hz))
    else:
        return 10**f2(np.log10(freq_hz))


def system_temp(freq_hz):
    """Return SKA1 Low system temperatures for a given frequency.

    Values provided by Eloy de Lera Acedo
    (email: eloy .at. mrao.cam.ac.uk)

    Args:
        freq_hz (float): Frequency, in Hz

    Returns:
        System temperature, in K
    """
    freqs, t_sys = np.array(raw_t_sys['freqs']), np.array(raw_t_sys['values'])
    f = interp1d(np.log10(freqs), np.log10(t_sys), kind='cubic')
    return 10**f(np.log10(freq_hz))


def evaluate_noise_rms(freq_hz, t_acc=5.0, bw_hz=100e3, eta=1.0,
                       num_antennas=256):
    """Evaluate the station Stokes-I RMS noise of a station, in Jy

    Args:
        freq_hz (float): Frequency in, Hz
        t_acc (float): Integration time, in seconds
        bw_hz (float): Bandwidth, in Hz
        eta (float): System efficiency factor
        num_antennas(int): Number of antennas in a station.

    Returns:
        Noise RMS in Jy
    """
    t_sys = system_temp(freq_hz)
    a_eff = element_effective_area(freq_hz) * num_antennas

    # Single receiver polarisation SEFD
    sefd = (2.0 * const.k_B.value * t_sys * eta) / a_eff
    sigma_pq = (sefd * 1e26) / (2.0 * bw_hz * t_acc)**0.5
    # Stokes-I noise is from two receiver polarisations -> scale by 1 / sqrt(2)
    sigma_pq /= 2**0.5
    return sigma_pq


def evaluate_scaled_noise_rms(freq_hz, num_times, bw_hz=100e3, eta=1.0,
                              obs_length_h=1000.0, num_antennas=256):
    """Evaluate the station Stokes-I noise RMS, in Jy for the given observation
       length sampled with the specified number of time samples.

    Args:
        freq_hz (float): Frequency in, Hz
        num_times (int): Number of time samples
        bw_hz (float): Bandwidth, in Hz
        eta (float): System efficiency factor
        obs_length_h (float): Target observation length, in hours
        num_antennas(int): Number of antennas in a station.

    Returns:
        Noise RMS in Jy
    """
    t_acc = (obs_length_h * 3600) / num_times
    t_sys = system_temp(freq_hz)
    a_eff = element_effective_area(freq_hz) * num_antennas

    # Single receiver polarisation SEFD
    sefd = (2.0 * const.k_B.value * t_sys * eta) / a_eff
    sigma_pq = (sefd * 1e26) / (2.0 * bw_hz * t_acc)**0.5
    # Stokes-I noise is from two receiver polarisations -> scale by 1 / sqrt(2)
    sigma_pq /= 2**0.5
    return sigma_pq


def plot_telescope(x, y, r_cut, station_d=35.0, filename=None):
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
    if filename:
        dir_name = os.path.dirname(filename)
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)
        fig.savefig(filename)
    else:
        plt.show()
    plt.close(fig)


def plot_uv_coords(uu, vv, plot_r=None, units='m', filename=None, r_lim=None):
    fig, ax = plt.subplots(figsize=(8, 8), nrows=1, ncols=1)
    ax.set_aspect('equal')
    alpha_ = 0.01 if not r_lim or r_lim > 500.0 else 0.05
    ax.plot(uu, vv, 'k.', ms=2.0, alpha=alpha_)
    ax.plot(-uu, -vv, 'k.', ms=2.0, alpha=alpha_)
    if plot_r is not None:
        for r in plot_r:
            c = plt.Circle((0, 0), r, fill=False, color='r', lw=2.0, alpha=1.0)
            ax.add_artist(c)
    ax.grid()
    ax.set_xlabel('uu (%s)' % units)
    ax.set_ylabel('vv (%s)' % units)
    if not r_lim:
        r = (uu**2 + vv**2)**0.5
        r_lim = r.max() * 1.05
    ax.set_xlim(-r_lim, r_lim)
    ax.set_ylim(-r_lim, r_lim)
    if filename:
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        if not filename.endswith('.png'):
            filename += '.png'
        fig.savefig(filename)
    else:
        plt.show()
    plt.close(fig)


def load_telescope(r_cut, station_d, lon, lat, plot_filename=None):
    # Load telescope model
    coords = np.loadtxt(join('v5.tm', 'layout.txt'))
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r = (x**2 + y**2)**0.5

    x = x[np.where(r < r_cut)]
    y = y[np.where(r < r_cut)]
    z = z[np.where(r < r_cut)]
    num_stations = x.shape[0]
    if plot_filename:
        plot_telescope(x, y, r_cut, station_d=station_d, filename=plot_filename)
    x, y, z = convert_enu_to_ecef(x, y, z, lon, lat)
    return x, y, z


def generate_uvw_coords_m(x, y, z, obs_length_h, mjd_mid, t_acc, ra, dec,
                          filename=None):
    num_times = int((obs_length_h * 3600.0) / t_acc)
    mjd_start = mjd_mid - ((obs_length_h / 2) / 24)
    num_baselines = x.shape[0] * (x.shape[0] - 1) // 2
    num_vis = num_baselines * num_times
    uu, vv, ww = np.zeros(num_vis), np.zeros(num_vis), np.zeros(num_vis)
    print('- Simulating uv coordinates for %i times over %.1f h'
          % (num_times, obs_length_h))
    t0 = time.time()
    for i in range(num_times):
        mjd = mjd_start + (i * t_acc + t_acc / 2) / 86400
        uu_, vv_, ww_ = evaluate_baseline_uvw(x, y, z, ra, dec, mjd)
        uu[i * num_baselines:(i + 1) * num_baselines] = uu_
        vv[i * num_baselines:(i + 1) * num_baselines] = vv_
        ww[i * num_baselines:(i + 1) * num_baselines] = ww_
    print('- Generated uv %i coordinates in %.2f s'
          % (uu.shape[0], time.time() - t0))
    if filename:
        dir_name = os.path.dirname(filename)
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)
        plot_uv_coords(uu, vv, filename=filename + '_all')
        plot_uv_coords(uu, vv, filename=filename + '_r%06.1fm' % 100.0,
                       r_lim=100.0)
    return uu, vv, ww, num_times


def write_fits_cube(cube, filename, ra, dec, mjd_start, freq_start, freq_inc,
                    fov, lambda_cut, weights, algorithm, bunit='Jy/beam'):
    # TODO(BM) convert units jy/beam -> K
    # TODO(BM) add fits header, WCS etc.
    # TODO(BM) write out one plane at a time? http://goo.gl/owain9
    hdu = fits.PrimaryHDU(cube)

    size = cube.shape[-1]
    header = hdu.header
    lm_max = sin(radians(fov) * 0.5)
    lm_inc = (lm_max * 2) / size
    cdelt = degrees(asin(lm_inc))
    crpix = size / 2 + 1

    time_utc = time.gmtime()
    header.set('DATE-OBS', '%04i-%02i-%02i' %
               (time_utc.tm_year, time_utc.tm_mon, time_utc.tm_mday))

    header.set('OBJECT', 'simulation of random noise')
    header.set('TELESCOP', 'SKA1 low v5')
    header.set('BUNIT', bunit)
    header.set('BMIN', cube.min(), bunit)
    header.set('BMAX', cube.max(), bunit)
    header.set('OBSRA', ra, 'RA')
    header.set('OBSDEC', dec, 'DEC')
    header.set('MJD-OBS', mjd_start, 'Start of observation')

    # Image axes
    header.set('CTYPE1', 'RA--SIN', 'Right Ascension')
    header.set('CUNIT1', 'deg')
    header.set('CRVAL1', ra, 'coordinate value at ref point')
    header.set('CRPIX1', crpix, 'pixel coordinate of the reference point')
    header.set('CDELT1', -cdelt, 'coordinate increment')

    header.set('CTYPE2', 'DEC--SIN', 'Declination')
    header.set('CUNIT2', 'deg')
    header.set('CRVAL2', dec, 'coordinate value at ref point')
    header.set('CRPIX2', crpix, 'pixel coordinate of the reference point')
    header.set('CDELT2', cdelt, 'coordinate increment')

    # Frequency axis
    header.set('CTYPE3', 'FREQ', 'Frequency')
    header.set('CUNIT3', 'Hz')
    header.set('CRVAL3', freq_start + freq_inc /2,
               'coordinate value at ref point')
    header.set('CRPIX3', 1, 'pixel coordinate of the reference point')
    header.set('CDELT3', freq_inc, 'coordinate increment')

    header.set('comment', 'Simulation of random noise for testing the EoR '
                          'pipeline')
    header.set('comment', 'Generated on %s' % time.asctime())
    header.set('comment', 'FOV = %f deg.' % fov)
    header.set('comment', 'weights = %s' % weights)
    header.set('comment', 'algorithm = %s' % algorithm)
    if lambda_cut:
        header.set('comment', 'lambda cut inner = %i' % lambda_cut[0])
        header.set('comment', 'lambda cut outer = %i' % lambda_cut[1])

    hdu_list = fits.HDUList([hdu])
    hdu_list.writeto(filename, clobber=True)


def gauss2d(xy, sx, sy, theta):
    theta = radians(theta)
    a = (cos(theta)**2 / (2 * sx**2) + sin(theta)**2 / (2 * sy**2))
    b = (-sin(2 * theta) / (4 * sx**2) + sin(2 * theta) / (4 * sy**2))
    c = (sin(theta)**2 / (2 * sx**2) + cos(theta)**2 / (2 * sy**2))
    x, y = xy
    inner = a * x**2
    inner -= 2 * b * x * y
    inner += c * y**2
    return np.exp(-inner)


def eval_beam_area(psf_cube, l_cut_outer, fov_deg, freqs, start_freq,
                   results_dir, weights, plot=False):
    fit_plot_dir = join(results_dir, 'psf_fit')
    if plot and not os.path.isdir(fit_plot_dir):
        os.makedirs(fit_plot_dir)

    num_channels = psf_cube.shape[0]
    im_size = psf_cube.shape[1]

    centre = im_size // 2
    lm_max = sin(radians(fov_deg) * 0.5)
    lm_inc = (2 * lm_max) / im_size
    hpbw_lm = sin(1 / l_cut_outer)
    sigma_lm = hpbw_lm / (2 * (2 * log(2))**0.5)

    sigma_x_arcmin = np.zeros(num_channels)
    sigma_y_arcmin = np.zeros(num_channels)
    theta_deg = np.zeros(num_channels)
    area_arcmin2 = np.zeros(num_channels)
    fit_rms = np.zeros(num_channels)

    for i in range(num_channels):
        image = psf_cube[i, :, :]
        # Crop to the ~first null
        vmin = 1e-5
        i0 = np.argmax(image[centre, centre:] < vmin)
        i0 = int(i0 * 1.25)
        c1 = centre - i0
        c2 = centre + i0 + 1
        crop = image[c1:c2, c1:c2]

        # Get grop image coordinates
        crop_size = crop.shape[0]
        crop_centre = crop_size // 2
        l = np.arange(-crop_centre, crop_centre + 1) * lm_inc
        lg, mg = np.meshgrid(-l, l)
        xg, yg = lg, mg
        xy = np.vstack((xg.flatten(), yg.flatten()))
        x, y, = xy

        guess = [sigma_lm, sigma_lm, 0.0]
        pred_params, uncert_cov = opt.curve_fit(gauss2d, xy, crop.flatten(),
                                               p0=guess)
        zpred = gauss2d(xy, *pred_params)
        sigma_x_lm, sigma_y_lm = pred_params[0], pred_params[1]
        sigma_x_arcmin[i] = degrees(asin(sigma_x_lm)) * 60
        sigma_y_arcmin[i] = degrees(asin(sigma_y_lm)) * 60
        theta_deg[i] = pred_params[2]
        area_arcmin2[i] = 2 * pi * sigma_x_arcmin[i] * sigma_y_arcmin[i]
        fit_rms[i] = (np.mean((crop.flatten() - zpred)**2))**0.5
        # Plot of fit
        if plot:
            fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=1, ncols=2)
            ax1.plot(-l, crop[crop_centre, :], c='b', lw=1.0, label='PSF')
            ax1.plot(-l, zpred.reshape(crop_size, crop_size)[crop_centre, :],
                     ':', c='r', label='Fit')
            ax1.set_xlim(-l[0], -l[-1])
            ax1.set_ylim(vmin, 1.5)
            ax1.set_title('x cut : %.3f MHz' % (freqs[i] / 1e6))
            ax2.set_xlabel('l direction cosine')
            ax2.set_ylabel('PSF amplitude')
            ax1.set_yscale('log')
            ax1.legend()
            ax2.plot(l, crop[:, crop_centre], c='b', lw=1.0, label='PSF')
            ax2.plot(-l, zpred.reshape(crop_size, crop_size)[:, crop_centre],
                     ':', c='r', label='Fit')
            ax2.set_xlim(-l[0], -l[-1])
            ax2.set_ylim(vmin, 1.5)
            ax2.set_title('y cut : %.3f MHz' % (freqs[i] / 1e6))
            ax2.set_xlabel('m direction cosine')
            ax2.set_ylabel('PSF amplitude')
            ax2.set_yscale('log')
            ax2.legend()
            fig.savefig(join(fit_plot_dir, 'fit_xy_cut_%05.1fMHz_c%03i_%s.png'
                             % (start_freq / 1e6, i, weights)))
            plt.close(fig)

    psf_fit_file = join(results_dir, 'psf_fit_%05.1fMHz_%s.npz'
                        % (start_freq / 1e6, weights))
    np.savez_compressed(psf_fit_file, sigma_x_arcmin=sigma_x_arcmin,
                        sigma_y_arcmin=sigma_y_arcmin, theta_deg=theta_deg,
                        area_arcmin2=area_arcmin2, fit_rms=fit_rms, freqs=freqs,
                        start_freq=start_freq)
    return area_arcmin2


def main():
    # Options
    # =========================================================================
    # Telescope model
    lon, lat = 116.63128900, -26.69702400  # degrees
    station_d = 35.0  # m
    telescope_radius_cut = 1000  # m (Remote stations in the telescope > this)
    eta = 1.0  # System efficiency

    # Observation settings
    az0, el0, date0 = 0.0, 90.0, '2016/7/14 10:30'  # Midpoint coordinates
    #start_freqs = 50e6 + np.array([0, 6, 12, 18]) * 8e6
    # start_freqs = [50e6]
    start_freqs = 50e6 + np.array([6, 12, 18]) * 8e6
    num_channels, freq_inc = 80, 100e3
    obs_length_h, noise_obs_length_h = 5, 1000
    t_acc = 60.0  # Time resolution of visibility coordinates, in seconds.

    # Image settings
    im_size = 512
    algorithm = 'w-projection'
    weights = 'uniform'
    lambda_cut = True
    w_planes = 64

    # Results directory
    results_dir = 'noise_cubes_%s_%s' % (weights, algorithm.lower())
    if algorithm.lower().startswith('w') and w_planes > 0:
        results_dir += '_%i' % w_planes
    # =========================================================================

    # Calculate observation equatorial coordinates and start MJD of the
    # specified target field
    obs = ephem.Observer()
    obs.lon, obs.lat, obs.elevation = radians(lon), radians(lat), 0.0
    obs.date = str(date0)
    ra, dec = obs.radec_of(radians(az0), radians(el0))
    ra, dec = degrees(ra), degrees(dec)
    mjd_mid = ephem.julian_date(b'2016/02/25 10:30') - 2400000.5
    mjd_start = mjd_mid - obs_length_h / (2 * 24.0)

    # Load telescope model and generate uvw coordinates, in metres
    coords_file = join(results_dir, 'uvw_m.npz')
    x, y, z = load_telescope(telescope_radius_cut, station_d, lon, lat,
                             join(results_dir, 'telescope.eps'))
    num_stations = x.shape[0]
    uu, vv, ww, num_times = generate_uvw_coords_m(x, y, z, obs_length_h,
                                                  mjd_mid, t_acc, ra, dec,
                                                  join(results_dir, 'uvw'))
    t0 = time.time()
    np.savez(coords_file, uu=uu, vv=vv, ww=ww, x=x, y=y, z=z,
             num_times=num_times)
    print('- No. stations = %i' % num_stations)
    print('- Saved uvw coordinates in %.2f s' % (time.time() - t0))

    # Evaluate radial uv co-ordinate range (needed for lambda cut)
    r_uvw = (uu**2 + vv**2)**0.5
    r_uvw_min = r_uvw.min()
    r_uvw_max = r_uvw.max()

    # Allocate image cubes.
    noise_cube = np.zeros((num_channels, im_size, im_size))
    psf_cube = np.zeros((num_channels, im_size, im_size))

    # Loop over cube start frequency.
    for i, start_freq in enumerate(start_freqs):
        print()
        print('- Freq = %.1f MHz (+%i x %0.1f kHz)' %
              (start_freq / 1e6, num_channels, freq_inc/1e3))

        # Array of cube frequency channels
        freqs = (start_freq + (np.arange(num_channels) * freq_inc + freq_inc / 2))

        # Calculate image FoV and lambda cuts
        fov_deg = 1.5 * degrees((const.c.value / start_freq) / station_d)
        print('- FoV = %.2f deg.' % fov_deg)

        # Generate file names of image cube and noise outputs.
        if lambda_cut:
            l_cut_inner = ceil(r_uvw_min / (const.c.value / freqs[-1]))
            l_cut_outer = floor(r_uvw_max / (const.c.value / freqs[0]))
            print('- lambda cut = %.0f, %.0f' % (l_cut_inner, l_cut_outer))
            suffix = ('%05.1fMHz_lcut_%04i_%04i_%s_%s' %
                      (start_freq / 1e6, l_cut_inner, l_cut_outer, weights,
                       algorithm.lower()))
        else:
            suffix = ('%05.1fMHz_%s_%s' % (start_freq / 1e6, weights,
                                           algorithm.lower()))
        if algorithm.lower().startswith('w') and w_planes > 0:
            suffix += '_%i' % w_planes
        l_cut = [l_cut_inner, l_cut_outer] if lambda_cut else None

        noise_cube_file = join(results_dir, 'noise_' + suffix + '.fits')
        noise_cube_file_k = join(results_dir, 'noise_' + suffix + '_K.fits')
        psf_cube_file = join(results_dir, 'psf_' + suffix + '.fits')
        noise_sigma_file = join(results_dir, 'noise_sigma_' + suffix + '.npz')

        t0 = time.time()
        # Option to load existing image cubes (used to skip directly to fitting
        # the psf and converting to K)
        if os.path.exists(noise_cube_file):
            noise_cube = fits.getdata(noise_cube_file)
            psf_cube = fits.getdata(psf_cube_file)
        # Simulate / image the noise cube (in Jy/beam) and PSF cube.
        else:
            sigma_pq, sigma_im = np.zeros(num_channels), np.zeros(num_channels)
            num_coords = np.zeros(num_channels, dtype='i8')
            progress_bar = ProgressBar(maxval=num_channels, widgets=[
                Bar(marker='='), Counter(format='%03i'),
                '/%03i ' % num_channels, Percentage(), ' ', Timer(), ' ',
                ETA()]).start()
            # Loop over frequencies in the cube
            for j, freq in enumerate(freqs):
                # Convert coordinates wavelength and apply lambda cut
                wavelength = const.c.value / freq
                uu_l = (uu / wavelength)
                vv_l = (vv / wavelength)
                ww_l = (ww / wavelength)
                if lambda_cut:
                    r_uv = (uu_l**2 + vv_l**2)**0.5
                    idx = np.where(np.logical_and(r_uv >= l_cut_inner,
                                                  r_uv <= l_cut_outer))
                    uu_l, vv_l, ww_l = uu_l[idx], vv_l[idx], ww_l[idx]
                num_coords[j] = uu_l.shape[0]

                # Evaluate visibility noise (and predicted image noise)
                sigma_pq[j] = evaluate_scaled_noise_rms(
                    freq, num_times, freq_inc, eta, noise_obs_length_h)
                sigma_im[j] = sigma_pq[j] / num_coords[j]**0.5

                # Make the noise image
                amp = (randn(num_coords[j]) + 1j * randn(num_coords[j])) * sigma_pq[j]
                noise_cube[j, :, :] = Imager.make_image(uu_l, vv_l, ww_l,
                                                        amp,
                                                        fov_deg, im_size,
                                                        weights, algorithm,
                                                        wprojplanes=w_planes)

                # Make the psf
                amp = np.ones(uu_l.shape[0], dtype='c16')
                psf_cube[j, :, :] = Imager.make_image(uu_l, vv_l, ww_l, amp,
                                                      fov_deg, im_size, weights,
                                                      algorithm,
                                                      wprojplanes=w_planes)
                progress_bar.update(j)
            progress_bar.finish()

            # Save noise values and cube images in Jy/beam
            np.savez(noise_sigma_file, freqs=freqs, sigma_pq=sigma_pq,
                     sigma_im=sigma_im, num_coords=num_coords)
            write_fits_cube(noise_cube, noise_cube_file, ra, dec, mjd_start,
                            start_freq, freq_inc, fov_deg,
                            l_cut, weights, algorithm)
            write_fits_cube(psf_cube, psf_cube_file, ra, dec, mjd_start,
                            start_freq, freq_inc, fov_deg,
                            l_cut, weights, algorithm)
        print('- Time taken = %.2f s (image cube)' % (time.time() - t0))

        # Fit the PSF with a gaussian to evaluate the beam area (in arcmin^2)
        t0 = time.time()
        r_uv_max_l = l_cut_outer if lambda_cut else \
            r_uvw_max / (const.c.value / start_freq)
        area_arcmin2 = eval_beam_area(psf_cube, r_uv_max_l, fov_deg, freqs,
                                      start_freq, results_dir, weights,
                                      plot=False)
        print('- Time taken = %.2f s (fit area)' % (time.time() - t0))

        # Convert cube from Jy/beam to K
        t0 = time.time()
        for j, freq in enumerate(freqs):
            image = noise_cube[j, :, :]
            area_sr = (area_arcmin2[j] / 3600.0) * (pi / 180)**2
            image /= area_sr
            image *= (1e-26 * const.c.value**2) / (2 * const.k_B.value * freq**2)
        write_fits_cube(noise_cube, noise_cube_file_k, ra, dec, mjd_start,
                        start_freq, freq_inc, fov_deg,
                        l_cut, weights, algorithm,
                        bunit='K')
        print('- Time taken = %.2f s (convert to K)' % (time.time() - t0))
        print()


if __name__ == '__main__':
    main()





