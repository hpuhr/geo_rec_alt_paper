#!/usr/bin/python3

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import util.common as common

def plot_adsb_geo_baro_offsets():

    geobaro_offsets = (common.cat021_df['geometric_height'] - common.cat021_df['mode_c_code']) * common.FT2M # diff ft to m

    plt.hist(geobaro_offsets, bins=common.num_bins, log=True, color='r')
    plt.xlabel('Geo-Baro Offset [m]')
    plt.ylabel('Count')
    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.savefig(common.output_folder + '/' + 'geo-baro_histo', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()

    heatmap, xedges, yedges = np.histogram2d(geobaro_offsets, common.cat021_df['mode_c_code'] / 100, bins=common.num_bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.clf()
    ax = plt.gca()
    im = plt.imshow(heatmap.T, interpolation='nearest', extent=extent, origin='lower', norm=LogNorm(), cmap='autumn')
    plt.xlabel('Geo-Baro Offset [m]')
    plt.ylabel('Barometric Altitude [FL]')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)

    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.savefig(common.output_folder + '/' + 'geo-baro_2d_histo', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()

def plot_cat021_geo_height():
    heatmap, xedges, yedges = np.histogram2d(
        common.cat021_df['geometric_height'] / 100.0, common.cat021_df['mode_c_code'] / 100.0, bins=common.num_bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.clf()
    ax = plt.gca()
    im = plt.imshow(heatmap.T, interpolation='nearest', extent=extent, origin='lower', norm=LogNorm())
    plt.xlabel('Geometric Altitude [FL]')
    plt.ylabel('Barometric Altitude [FL]')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)

    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.savefig(common.output_folder + '/' + 'cat021_geoheight_2d_histo', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()
