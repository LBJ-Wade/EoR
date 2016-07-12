# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
from math import cos, sin
from astropy.io import fits
from astropy import constants as const
from os.path import join
import os
from math import sin, radians, degrees, log, pi, asin
import seaborn
seaborn.set_style('ticks')


def gauss_param(theta, sigma_x, sigma_y):
    a = (cos(theta) ** 2 / (2 * sigma_x ** 2) +
         sin(theta) ** 2 / (2 * sigma_y ** 2))
    b = (-sin(2 * theta) / (4 * sigma_x ** 2) +
         sin(2 * theta) / (4 * sigma_y ** 2))
    c = (sin(theta) ** 2 / (2 * sigma_x ** 2) +
         cos(theta) ** 2 / (2 * sigma_y ** 2))
    return a, b, c


def main():
    x0, y0 = 0.0, 0.0
    amp = 1
    theta = np.radians(45)
    sigma_x = 0.2
    sigma_y = 1.0
    a, b, c = gauss_param(45, 0.2, 1.0)
    true_params = [amp, x0, y0, a, b, c]
    print(true_params)
    xy, zobs = generate_example_data(1000, true_params)
    x, y = xy

    # fig, ax = plt.subplots()
    # ax.scatter(x, y, c=zobs, s=50)
    # plt.show()

    i = zobs.argmax()
    pa, pb, pc = gauss_param(0, 0.4, 0.7)
    guess = [1, 0, 0, pa, pb, pc]
    pred_params, uncert_cov = opt.curve_fit(gauss2d, xy, zobs, p0=guess)

    zpred = gauss2d(xy, *pred_params)
    np.set_printoptions(precision=2)
    print('True parameters : ', true_params)
    print('Guess params    :', guess)
    print('Predicted params:', pred_params)
    print('Residual, RMS(obs - pred):', np.sqrt(np.mean((zobs - zpred)**2)))

    plot(xy, zobs, pred_params)
    plt.show()


def gauss2d(xy, amp, x0, y0, a, b, c):
    x, y = xy
    inner = a * (x - x0)**2
    inner -= 2 * b * (x - x0) * (y - y0)
    inner += c * (y - y0)**2
    return amp * np.exp(-inner)


def gauss2d_1(xy, a, b, c):
    x, y = xy
    x0, y0 = 0, 0
    inner = a * (x - x0)**2
    inner -= 2 * b * (x - x0) * (y - y0)
    inner += c * (y - y0)**2
    return np.exp(-inner)


def gauss2d_2(xy, sx, sy, theta):
    a, b, c = gauss_param(theta, sx, sy)
    return gauss2d_1(xy, a, b, c)


def generate_example_data(num, params):
    # xy = np.random.random((2, num)) * 4 - 2
    c2 = 10
    x = np.arange(-c2, c2 + 1) * 0.2
    xg, yg = np.meshgrid(x, x)
    xy = np.vstack((xg.flatten(), yg.flatten()))
    x, y = xy
    zobs = gauss2d(xy, *params)
    return xy, zobs


def plot(xy, zobs, pred_params):
    x, y = xy
    yi, xi = np.mgrid[-2:2:100j, -2:2:100j]
    xyi = np.vstack([xi.ravel(), yi.ravel()])

    zpred = gauss2d(xyi, *pred_params)
    zpred.shape = xi.shape

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=zobs, s=50, vmin=zpred.min(), vmax=zpred.max())
    im = ax.imshow(zpred, extent=[xi.min(), xi.max(), yi.max(), yi.min()],
                   aspect='auto', interpolation='nearest')
    fig.colorbar(im)
    cs = ax.contour(xi, yi, zpred, [0.1, 0.5], colors='w', linewidths=[1.0, 2.0],
               linestyles='-')
    ax.clabel(cs, inline=1, fontsize=10)
    ax.invert_yaxis()
    return fig


def main2():
    data = fits.getdata(join('results', 'noise_050.0-050.3MHz_100kHz_5h_60s_natural',
                             'psf_no_l_cut_natural.fits'))
    data = data[0, :, :]
    c = data.shape[0]//2
    i0 = np.where(data[c, c:] < 0.0)[0][0]
    i0 = int(i0 * 1.5)
    # i0 = 2
    data2 = data[c - i0:c + i0 + 1, c - i0:c + i0 + 1]
    freq = 50e6 + (np.arange(200) * 100e3 + 50e3)

    # Get lm cell size for the entire image
    fov = 10.0  # degrees
    size = data.shape[0]
    centre = size // 2
    lm_max = sin(radians(fov) * 0.5)
    lm_inc = (2 * lm_max) / size
    # extent for the entire image
    extent = np.array([centre + 0.5, -centre + 0.5,
                       -centre - 0.5, centre - 0.5])
    extent *= lm_inc

    # Get the extent for the inner cut reguion, 'data2'
    size2 = data2.shape[0]
    centre2 = size2 // 2
    l = np.arange(-centre2, centre2 + 1) * lm_inc
    lg, mg = np.meshgrid(-l, l)
    extent2 = np.array([centre2 + 0.5, -centre2 - 0.5,
                        -centre2 - 0.5, centre2 + 0.5])
    extent2 *= lm_inc

    hpbw = degrees(1 / 230)
    hpbw_lm = sin(1 / 230)
    print(hpbw, hpbw_lm)

    xg, yg = lg, mg
    x, y = -l, l
    xy = np.vstack((xg.flatten(), yg.flatten()))
    x, y, = xy

    # print(x.shape, y.shape, data2.flatten().shape)
    # fig, ax = plt.subplots()
    # ax.scatter(x, y, c=data2.flatten(), s=20, alpha=0.8, lw=0.5)
    # ax.set_aspect('equal')
    # plt.show()

    pa, pb, pc = gauss_param(theta=0, sigma_x=hpbw_lm/2, sigma_y=hpbw_lm/2)
    guess = [1, 0, 0, pa, pb, pc]
    pred_params, uncert_cov = opt.curve_fit(gauss2d, xy, data2.flatten(),
                                            p0=guess)
    zpred = gauss2d(xy, *pred_params)
    np.set_printoptions(precision=2)
    print('guess    :', guess[3:])
    print('predicted:', pred_params[3:])
    print('Residual, RMS(obs - pred):',
          np.sqrt(np.mean((data2.flatten() - zpred)**2)))

    fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=1, ncols=2)
    # axis1
    im = ax1.imshow(data2, interpolation='nearest', origin='lower',
                    extent=extent2, alpha=1.0, cmap='inferno',
                    vmin = -0.05, vmax = 0.6)
    c = plt.Circle((0, 0), hpbw_lm / 2, fill=False, color='w', lw=2.0)
    ax1.add_artist(c)
    cs = ax1.contour(xg, yg, data2, [0.5], colors='c', lw=2.0)
    ax1.clabel(cs, inline=1, fontsize=10)
    ax1.set_title('psf')
    # axis2
    im = ax2.imshow(zpred.reshape(size2, size2), extent=extent2,
                    interpolation='nearest', origin='lower',
                    alpha=1.0, cmap='inferno', vmin=-0.05, vmax=0.6)
    cs = ax2.contour(xg, yg, zpred.reshape(size2, size2), [0.5],
                     colors='c', lw=2.0)
    ax2.clabel(cs, inline=1, fontsize=10)
    c = plt.Circle((0, 0), hpbw_lm / 2, fill=False, color='w', lw=2.0)
    ax2.add_artist(c)
    ax2.set_xlim(extent2[0], extent2[1])
    ax2.set_ylim(extent2[2], extent2[3])
    ax2.set_aspect('equal')
    ax2.set_title('gaussian fit')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.2, 0.03, 0.6])
    fig.colorbar(im, cax=cbar_ax)

    plt.show()


def main3():
    # TODO(BM) take the log and fit a parabola instead of the gaussian fit used
    # here https://en.wikipedia.org/wiki/Gaussian_function
    # (see profile estimation)

    # weight = 'uniform'
    weight = 'natural'
    freq0 = 200e6

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
    cube_file = join('results',
                     'noise_%05.1f-%05.1fMHz_100kHz_5h_60s_%s' %
                     (freq0/1e6, freq1/1e6, weight),
                     'psf_l_cut_%i_%i_%s.fits' %
                     (b_min_lambda, b_max_lambda, weight))
    cube = fits.getdata(cube_file)

    out_dir = os.path.dirname(cube_file)

    # TODO(BM) get frequency info from FITS header.
    freq = freq0 + (np.arange(200) * 100e3 + 50e3)

    # Get lm cell size for the entire image
    # TODO(BM) get from FITS header
    fov = 10.0  # degrees
    size = cube.shape[1]
    centre = size // 2
    lm_max = sin(radians(fov) * 0.5)
    lm_inc = (2 * lm_max) / size
    # extent for the entire image
    extent = np.array([centre + 0.5, -centre + 0.5,
                       -centre - 0.5, centre - 0.5])
    extent *= lm_inc
    num_channels = cube.shape[0]

    # TODO(BM) get from fits header
    hpbw = degrees(1 / b_max_lambda)
    hpbw_lm = sin(1 / b_max_lambda)
    sigma = hpbw / (2*(2*log(2))**0.5)
    sigma_lm = hpbw_lm / (2*(2*log(2))**0.5)

    plot_1d = True
    plot_2d = False
    plot_dir = join(out_dir, 'psf_fit')
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)

    area = np.zeros(num_channels)
    sigma_x_arcmin = np.zeros(num_channels)
    sigma_y_arcmin = np.zeros(num_channels)
    theta_deg = np.zeros(num_channels)
    fit_rms = np.zeros(num_channels)
    volume = np.zeros(num_channels)

    for i in range(num_channels):
        # Image plane
        image = cube[i, :, :]

        # Crop to the first null
        c = image.shape[0] // 2
        i0 = np.where(image[c, c:] < 0.0)[0][0]
        i0 = int(i0 * 1.5)
        image_crop = image[c - i0:c + i0 + 1, c - i0:c + i0 + 1]

        # Get the extent for the inner cut region, 'data2'
        size2 = image_crop.shape[0]
        centre2 = size2 // 2
        l = np.arange(-centre2, centre2 + 1) * lm_inc
        lg, mg = np.meshgrid(-l, l)
        extent2 = np.array([centre2 + 0.5, -centre2 - 0.5,
                            -centre2 - 0.5, centre2 + 0.5])
        extent2 *= lm_inc

        xg, yg = lg, mg
        xy = np.vstack((xg.flatten(), yg.flatten()))
        x, y, = xy

        sigma_guess_lm = hpbw_lm / (2*(2*log(2))**0.5)
        guess = [sigma_guess_lm, sigma_guess_lm, 0.0]
        pred_params, uncert_cov = opt.curve_fit(gauss2d_2, xy,
                                                image_crop.flatten(),
                                                p0=guess)
        zpred = gauss2d_2(xy, *pred_params)
        np.set_printoptions(precision=5)

        sigma_x_lm = pred_params[0]
        sigma_y_lm = pred_params[1]
        sigma_x_arcmin[i] = degrees(asin(sigma_x_lm)) * 60.0
        sigma_y_arcmin[i] = degrees(asin(sigma_y_lm)) * 60.0
        theta_deg[i] = pred_params[2]
        area[i] = 2 * pi * sigma_x_arcmin[i] * sigma_y_arcmin[i]
        fit_rms[i] = np.sqrt(np.mean((image_crop.flatten() - zpred) ** 2))
        print('*' * 50)
        print('%i %.2f MHz' % (i, freq[i]/1e6))
        print('guess  : sx:%.2e, sy:%.2e, theta=%.1f' % tuple(guess))
        print('fitted : sx:%.2e, sy:%.2e, theta=%.1f' % tuple(pred_params))
        print('guess  : %f' % (degrees(asin(sigma_guess_lm)) * 60.0))
        print('sigma  : %f, %f' % (sigma_x_arcmin[i], sigma_y_arcmin[i]))
        # Cramér–Rao bound?
        # see section on Gaussian profile estimation @ https://goo.gl/X7CBDQ
        print('Residual, RMS(obs - pred):', fit_rms[i])
        print('beam area: %f' % area[i])
        print('guess area: %f ' %
              (2*pi*(degrees(asin(sigma_guess_lm)) * 60.0)**2))

        vmin_ = -0.1
        vmax_ = 1.0
        if plot_2d:
            fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=1, ncols=2)
            # axis1
            im = ax1.imshow(image_crop, interpolation='nearest', origin='lower',
                            extent=extent2, alpha=1.0, cmap='inferno',
                            vmin=vmin_, vmax=vmax_)
            c = plt.Circle((0, 0), hpbw_lm / 2, fill=False, color='w', lw=2.0)
            ax1.add_artist(c)
            cs = ax1.contour(xg, yg, image_crop, [0.5], colors='c', lw=2.0)
            ax1.clabel(cs, inline=1, fontsize=10)
            ax1.set_title('psf')
            # axis2
            im = ax2.imshow(zpred.reshape(size2, size2), extent=extent2,
                            interpolation='nearest', origin='lower',
                            alpha=1.0, cmap='inferno', vmin=vmin_, vmax=vmax_)
            cs = ax2.contour(xg, yg, zpred.reshape(size2, size2), [0.5],
                             colors='c', lw=2.0)
            ax2.clabel(cs, inline=1, fontsize=10)
            c = plt.Circle((0, 0), hpbw_lm / 2, fill=False, color='w', lw=2.0)
            ax2.add_artist(c)
            ax2.set_xlim(extent2[0], extent2[1])
            ax2.set_ylim(extent2[2], extent2[3])
            ax2.set_aspect('equal')
            ax2.set_title('gaussian fit')
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.2, 0.03, 0.6])
            fig.colorbar(im, cax=cbar_ax)
            fig.savefig(join(plot_dir, '2d_fit_%03i.png' % i))
            plt.close(fig)

        if plot_1d:
            # Plot x and y 1d cuts
            fig, (ax1, ax2) = plt.subplots(figsize=(10, 6), nrows=1, ncols=2)
            # axis1 - x cut
            ax1.plot(-l, image_crop[centre2, :], label='psf', c='b')
            ax1.plot([hpbw_lm / 2, hpbw_lm / 2], [vmin_, vmax_], c='k',
                     label='guess fhwm')
            ax1.plot([-hpbw_lm / 2, -hpbw_lm / 2], [vmin_, vmax_], c='k')
            ax1.plot([-(2 * log(2))**0.5 * sigma_x_lm,
                      -(2 * log(2))**0.5 * sigma_x_lm],
                     [vmin_, vmax_], '--', c='g', label='fit fwhm')
            ax1.plot([(2 * log(2))**0.5 * sigma_x_lm,
                      (2 * log(2))**0.5 * sigma_x_lm],
                     [vmin_, vmax_], '--', c='g')
            ax1.plot(-l, zpred.reshape(size2, size2)[centre2, :],
                     label='gauss fit', c='r')
            ax1.plot(ax1.get_xlim(), [0.5, 0.5], '--', c='0.5')
            ax1.grid()
            ax1.set_title('x cut - %.3f MHz' % (freq[i]/1e6))
            ax1.set_xlim(-l[0], -l[-1])
            ax1.legend()
            # axis2 - y cut
            ax2.plot(l, image_crop[:, centre2], label='psf', c='b')
            ax2.plot([hpbw_lm / 2, hpbw_lm / 2], [vmin_, vmax_], c='k',
                     label='guess fwhm')
            ax2.plot([-hpbw_lm / 2, -hpbw_lm / 2], [vmin_, vmax_], c='k')
            ax2.plot(l, zpred.reshape(size2, size2)[:, centre2],
                     label='gauss fit', c='r')
            ax2.plot([-(2 * log(2))**0.5 * sigma_y_lm,
                      -(2 * log(2))**0.5 * sigma_y_lm],
                     [vmin_, vmax_], '--', c='g', label='fit fwhm')
            ax2.plot([(2 * log(2))**0.5 * sigma_y_lm,
                      (2 * log(2))**0.5 * sigma_y_lm],
                     [vmin_, vmax_], '--', c='g')
            ax2.plot(ax2.get_xlim(), [0.5, 0.5], '--', c='0.5')
            ax2.set_title('y cut - %.3f MHz' % (freq[i] / 1e6))
            ax2.set_xlim(l[0], l[-1])
            ax2.legend()
            ax2.grid()
            fig.savefig(join(plot_dir, '1d_fit_%03i.png' % i))
            plt.close(fig)

    psf_fit_data = np.vstack((sigma_x_arcmin, sigma_y_arcmin, theta_deg,
                              fit_rms, area))
    header_ = 'sigma_x (arcmin), sigma_y (arcmin), theta (deg), ' \
              'fit rms, area (arcmin^2)'
    np.savetxt(join(plot_dir, 'fit_data.txt'), np.transpose(psf_fit_data),
               header=header_)

if __name__ == '__main__':
    # main()
    # main2()
    main3()
