#!/usr/bin/python3

import math
import pygeodesy
import multiprocessing
import pandas as pd
import matplotlib.pyplot as plt
import random
import statistics

import util.common as common
from reconstructor_cell import ReconstructorCell

def create():
    # add height offsets for calculation
    print("calculating adsb geometric height offsets")
    ginterpolator = pygeodesy.GeoidPGM(common.egm96_filename)
    common.cat021_df['geometric_height_correction_m'] = ginterpolator.height(
        common.cat021_df['latitude'], common.cat021_df['longitude'])

    print("inializing adsb reconstructor cells")
    common.cat021_df["latitude_bins"] = pd.cut(
        common.cat021_df["latitude"], common.latitude_bins, labels=range(len(common.latitude_bins) - 1))
    common.cat021_df["longitude_bins"] = pd.cut(
        common.cat021_df["longitude"], common.longitude_bins, labels=range(len(common.longitude_bins) - 1))
    common.cat021_df["time_of_day_bins"] = pd.cut(
        common.cat021_df["time_of_day"], common.time_of_day_bins, labels=range(len(common.time_of_day_bins) - 1))

    # print(cat021_df["longitude_bins"])

    adsb_reconst_cells = [[[ReconstructorCell(tod_bin, common.time_of_day_bins,
                                              lat_bin, common.latitude_bins,
                                              lon_bin, common.longitude_bins)
                            for lon_bin in range(len(common.longitude_bins))]
                           for lat_bin in range(len(common.latitude_bins))]
                          for tod_bin in range(len(common.time_of_day_bins))]

    print("filling adsb reconstructor cells")
    for cnt in range(common.cat021_df.shape[0]):

        if cnt % 200000 == 0:
            print('creating indexes {}'.format(cnt))

        # data bins
        lat_bin = common.cat021_df["latitude_bins"].iloc[cnt]
        lon_bin = common.cat021_df["longitude_bins"].iloc[cnt]
        tod_bin = common.cat021_df["time_of_day_bins"].iloc[cnt]

        if math.isnan(tod_bin):
            print('tod_bin {} tod {}'.format(tod_bin, common.cat021_df["time_of_day"].iloc[cnt]))
            continue

        if math.isnan(lat_bin):
            print('lat_bin {} lat {}'.format(lat_bin, common.cat021_df["latitude"].iloc[cnt]))
            continue

        if math.isnan(lon_bin):
            print('lon_bin {} lon {}'.format(lon_bin, common.cat021_df["longitude"].iloc[cnt]))
            continue

        adsb_reconst_cells[tod_bin][lat_bin][lon_bin].indexes.append(cnt)

    print('processing cells with {} tr'.format(common.cat021_df.shape[0]))

    num_cells = 0

    max_cell_size = None

    # postprocess cells
    for time_cells in adsb_reconst_cells:
        for lat_index in range(len(time_cells)):
            # for lat_cells in time_cells:
            lat_cells = time_cells[lat_index]

            pool = multiprocessing.Pool()  # initialise your pool
            finalized_cells = pool.map(ReconstructorCell.finalize, lat_cells)
            pool.close()  # shut down the pool
            pool.join()

            for cell in finalized_cells:  # type: ReconstructorCell

                num_cells += 1

                if len(cell.indexes) > 0:
                    #cells.append(cell)

                    if max_cell_size is None or len(cell.indexes) > max_cell_size:
                        max_cell_size = len(cell.indexes)

            # finalized_cells
            time_cells[lat_index] = finalized_cells

    print('cells {} avg size {} max size {}'.format(num_cells, common.cat021_df.shape[0] / num_cells, max_cell_size))

    return adsb_reconst_cells

def calc_reconst_offsets(ratio):
    num_misses = 0
    offsets_equal = []
    offsets_reconst = []

    data_bins_len = len(common.cat021_df["time_of_day_bins"])

    num_tests = int(data_bins_len * ratio)

    test_indexes = random.sample(list(range(data_bins_len)), num_tests)

    print('testing {} samples from {} target reports'.format(len(test_indexes), data_bins_len))

    for test_index in test_indexes:

        tod_bin = common.cat021_df["time_of_day_bins"].iloc[test_index]
        lat_bin = common.cat021_df["latitude_bins"].iloc[test_index]
        lon_bin = common.cat021_df["longitude_bins"].iloc[test_index]

        # print('tod_bin {} lat_bin {} lon_bin {}'.format(tod_bin, lat_bin, lon_bin))

        mc_ft = common.cat021_df['mode_c_code'].iloc[test_index]
        geo_ft = common.cat021_df['geometric_height'].iloc[test_index]

        assert not math.isnan(mc_ft)
        assert not math.isnan(geo_ft)

        cell = common.adsb_reconst_cells[tod_bin][lat_bin][lon_bin] # type: ReconstructorCell

        if not cell.param_estimated:
            num_misses += 1
            continue

        alt_cor = cell.calc_alt_corr(mc_ft)

        offsets_equal.append(abs(geo_ft - mc_ft) * common.FT2M)
        offsets_reconst.append(abs(geo_ft - alt_cor) * common.FT2M)

    total_tests = num_misses + len(offsets_reconst)

    print('cell sizes tod {} lat/lon {}'.format(common.tod_step, common.geo_step))
    print('test performed {} misses {} perc {}'.format(total_tests, num_misses, 100 * num_misses / total_tests))
    print('test equal errors avg {} std.dev. {}'.format(
        statistics.mean(offsets_equal), statistics.stdev(offsets_equal)))
    print('test reconst errors avg {} std.dev. {}'.format(
        statistics.mean(offsets_reconst), statistics.stdev(offsets_reconst)))

    return num_misses, offsets_equal, offsets_reconst



def plot_reconst_errors(offsets_equal, offsets_reconst):

    xmin = min(min(offsets_equal), min(offsets_reconst))
    xmax = max(max(offsets_equal), max(offsets_reconst))

    x_range = [xmin, xmax]
    common.num_bins = 50

    # all
    plt.hist(offsets_equal, range=x_range, bins=common.num_bins, log=True, alpha=0.5, label='Geo-Baro', color='r')
    plt.hist(offsets_reconst, range=x_range, bins=common.num_bins, log=True, alpha=0.5, label='Geo-Reconst', color='b')
    plt.xlabel('Altitude Error [m]')
    plt.ylabel("Count")
    plt.legend(loc='upper right')
    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.savefig(common.output_folder + '/' + 'reconst_errors_histo.png', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()
