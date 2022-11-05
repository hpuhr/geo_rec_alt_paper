#!/usr/bin/python3

import sys
import sqlite3
import time
import numpy as np

import util.common as common
import points_of_interest
import load_data
import adsb_cells
import adsb_plots
import radar_plots
import tracker_chains
import radar_analysis


def main(argv):
    start_time = time.time()

    common.latitude_bins = np.arange(common.latitude_min, common.latitude_max + common.geo_step, common.geo_step)
    common.longitude_bins = np.arange(common.longitude_min, common.longitude_max + common.geo_step, common.geo_step)
    common.time_of_day_bins = np.arange(0, 24 * 60 * 60 + common.tod_step, common.tod_step)

    common.create_path_if_required(common.output_folder)

    # Create a SQL connection to our SQLite database
    con = sqlite3.connect(common.db_filename)

    common.cat021_df = load_data.load_cat021(con)
    common.cat048_df = load_data.load_cat048(con)
    common.cat062_df = load_data.load_cat062(con)

    common.adsb_reconst_cells = adsb_cells.create()

    adsb_plots.plot_adsb_geo_baro_offsets()

    num_misses, offsets_equal, offsets_reconst = adsb_cells.calc_reconst_offsets(0.05)

    adsb_cells.plot_reconst_errors(offsets_equal, offsets_reconst)

    points_of_interest.calculate()

    radar_plots.plot_radar_3d_height()
    adsb_plots.plot_cat021_geo_height()

    common.trk_chains = tracker_chains.create()

    radar_analysis.test_coordinates()

    radar_analysis.calc_original_offsets()

    print('done after {}'.format(common.time_str_from_seconds(time.time() - start_time)))

    # Be sure to close the connection
    con.close()


if __name__ == "__main__":
    main(sys.argv[1:])

