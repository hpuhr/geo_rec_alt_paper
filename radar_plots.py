#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

import util.common as common

def plot_radar_3d_height():
    heatmap, xedges, yedges = np.histogram2d(
        common.cat048_df['height_3d'] / 100.0, common.cat048_df['mode_c_code'] / 100.0, bins=common.num_bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.clf()
    ax = plt.gca()
    im = plt.imshow(heatmap.T, interpolation='nearest', extent=extent, origin='lower', norm=LogNorm())
    plt.xlabel('3D Height [FL]')
    plt.ylabel('Barometric Altitude [FL]')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)

    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.savefig(common.output_folder + '/' + 'radar_3dheight_2d_histo', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()

