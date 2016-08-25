# -*- coding: utf-8 -*-
"""
http://balbuceosastropy.blogspot.co.uk/2013/09/the-mollweide-projection.html
"""
from __future__ import absolute_import, division, print_function
from mpl_toolkits.basemap import Basemap
import ephem
import matplotlib.pyplot as plt
import numpy as np


def plot_mwd(RA, Dec, org=0):
    x = np.remainder(RA+360-org,360) # shift RA values
    ind = x > 180
    x[ind] -= 360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection='mollweide', axisbg ='None')
    ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)

if __name__ == '__main__':
    lon, lat = 116.63128900, -26.69702400
    ra, dec = 68.698903779331502, -26.568851215532160

    fig = plt.figure(figsize=(10, 6))
    m = Basemap(projection='moll', lon_0=lon, resolution='c')
    m.drawcoastlines(linewidth=1.0)
    m.scatter(lon, lat, c='r', latlon=True)
    m.drawparallels(np.arange(-90, 120, 30), color='b')
    m.drawmeridians(np.arange(0, 420, 60), color='blue')
    #m.fillcontinents(color='w', lake_color='b')
    m.drawmapboundary(linewidth=1.0, fill_color='None')
    plt.show()

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='mollweide')
    ax.grid()
    ax.plot([ra], [dec], 'rx', ms=10, mew=2)
    plt.show()
    plot_mwd(np.array([ra]), np.array([dec]), org=116)
    plt.show()

    # lon_array = np.arange(0, 360)
    # lat = 0.
    # eq_array = np.zeros((360, 2))
    # for lon in lon_array:
    #     ga = ephem.Galactic(np.radians(lon), np.radians(lat))
    #     eq = ephem.Equatorial(ga)
    #     eq_array[lon] = np.degrees(eq.get())
    # RA = eq_array[:, 0]
    # Dec = eq_array[:, 1]
    # plot_mwd(RA, Dec, 116)
    # plt.show()
