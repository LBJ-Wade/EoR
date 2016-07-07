# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
from math import cos, sin
from astropy.io import fits
from astropy import constants as const
from os.path import join
from math import sin, radians, degrees
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


def generate_example_data(num, params):
    xy = np.random.random((2, num)) * 4 - 2
    # c2 = 10
    # x = np.arange(-c2, c2 + 1) * 0.2
    # xg, yg = np.meshgrid(x, x)
    # xy = np.vstack((xg.flatten(), yg.flatten()))
    # x, y = xy
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
    data = fits.getdata(join('results', 'psf_50-70MHz_100kHz_5h_60s_uniform',
                             'psf_l_cut_10_230_uniform.fits'))
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
                    vmin=-0.1, vmax=1.0)
    c = plt.Circle((0, 0), hpbw_lm / 2, fill=False, color='w', lw=2.0)
    ax1.add_artist(c)
    cs = ax1.contour(xg, yg, data2, [0.5], colors='c', lw=2.0)
    ax1.clabel(cs, inline=1, fontsize=10)
    ax1.set_title('psf')
    # axis2
    im = ax2.imshow(zpred.reshape(size2, size2), extent=extent2,
                    interpolation='nearest', origin='lower',
                    alpha=1.0, cmap='inferno', vmin=-0.1, vmax=1.0)
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


if __name__ == '__main__':
    # main()
    main2()
