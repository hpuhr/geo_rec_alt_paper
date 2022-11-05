#!/usr/bin/python3

import math
import matplotlib.pyplot as plt
import statistics

import numpy as np

import util.common as common

from numpy import linspace, random
from scipy.optimize import leastsq, least_squares

g0 = 9.81
R = 287

def residual(variables, ps, hs_data):
    """Model a decaying sine wave and subtract data."""
    P0 = variables[0]
    T0 = variables[1]
    alpha = variables[2]

    hs_model = T0 / alpha * (1 - (ps / P0) ** ((alpha * R) / g0))

    return hs_data-hs_model # / eps_data

def residual2(variables, ps, hs_data):
    """Model a decaying sine wave and subtract data."""
    P0 = variables[0]
    T0 = variables[1]
    alpha = 0.0065

    hs_model = T0 / alpha * (1 - (ps / P0) ** ((alpha * R) / g0))

    return hs_data-hs_model # / eps_data

class ReconstructorCell:

    def __init__(self, tod_bin, tod_bins, lat_bin, lat_bins, lon_bin, lon_bins):
        self.tod_bin = tod_bin
        self.tod_bins = tod_bins
        self.lat_bin = lat_bin
        self.lat_bins = lat_bins
        self.lon_bin = lon_bin
        self.lon_bins = lon_bins

        self.lat_min = None
        self.lat_max = None
        self.lon_min = None
        self.lon_max = None

        self.indexes = []

        self.finalized = False

        self.mcs_fts = []
        self.alt_geos = []

        self.calc_alt_geos = [] # corrected
        self.calc_ps = [] # from baro alt using std formula
        self.calc_p0s = [] # only used for scatterplot

        self.geo_alt_correction_avg_m = None

        self.param_estimated = False
        self.est_p0 = None
        self.est_t0 = None
        self.est_alpha = None

    def finalize(self):

        assert not self.finalized

        geo_alt_correction_avgs_m =[]

        for index in self.indexes:
            # print('mc {} geo {}'.format(v2_data['mode_c_fl'].iloc[index], v2_data['alt_geo_fl'].iloc[index]))

            lat = common.cat021_df['latitude'].iloc[index]
            lon = common.cat021_df['longitude'].iloc[index]

            if self.lat_min is None:
                self.lat_min = lat
                self.lat_max = lat

                self.lon_min = lon
                self.lon_max = lon
            else:
                self.lat_min = min(self.lat_min, lat)
                self.lat_max = max(self.lat_max, lat)

                self.lon_min = min(self.lon_min, lon)
                self.lon_max = max(self.lon_max, lon)

            mc_ft = common.cat021_df['mode_c_code'].iloc[index]
            alt_geo_ft = common.cat021_df['geometric_height'].iloc[index]
            geo_alt_correction_m = common.cat021_df['geometric_height_correction_m'].iloc[index]

            assert not math.isnan(mc_ft)
            assert not math.isnan(alt_geo_ft)

            # calculate pressure at mc for std formula wrt MSL (geoid)
            calc_p = common.calc_p_from_baro_std(mc_ft * common.FT2M)
            # calculate p0 that should have been used wrt MSL (geoid), gpsalt - geooidalt = msl height
            calc_p0 = common.calc_p0_from_p_geo(calc_p, alt_geo_ft * common.FT2M - geo_alt_correction_m)

            self.mcs_fts.append(mc_ft)
            self.alt_geos.append(alt_geo_ft)
            geo_alt_correction_avgs_m.append(geo_alt_correction_m)
            self.calc_alt_geos.append(alt_geo_ft * common.FT2M - geo_alt_correction_m) #
            self.calc_ps.append(calc_p)
            self.calc_p0s.append(calc_p0)

        self.finalized = True

        if len(self.indexes):
            self.geo_alt_correction_avg_m = statistics.mean(geo_alt_correction_avgs_m)

            self.do_param_reconstruction()

        return self

    def do_param_reconstruction(self):

        assert self.finalized

        if len(self.indexes) < 200:
            return

        # P0, T0, alpha
        variables = [1013.25, 273.15 + 15, 0.0065] #

        out = least_squares(residual, variables, args=(self.calc_ps, self.calc_alt_geos), method='lm')
        est_p0, est_t0, est_alpha = out.x #

        if est_p0 is not None and not math.isnan(est_p0) \
                and est_t0 is not None and not math.isnan(est_t0):

            self.est_p0 = est_p0
            self.est_t0 = est_t0
            self.est_alpha = est_alpha

            self.param_estimated = True

        return
        #return 'estimated P0 {:.2f} T0 {:.2f} alpha {:.5f}'.format(self.est_p0, self.est_t0 - 273.15, self.est_alpha)

    def write_calc(self, data_label_str, out_folder):

        assert self.finalized

        common.create_path_if_required(out_folder)

        file1 = open(out_folder + 'altitudes_cell_tod_{}_lat_{}_lon_{}.txt'.format(
            common.time_str_from_seconds(self.tod_bins[self.tod_bin]),
            self.lat_bins[self.lat_bin], self.lon_bins[self.lon_bin]), "w")

        file1.write('tod;mode_c_code;baro_m;geo_ft;geo_rec_ft;p\n')

        for index in self.indexes:
            mc_ft = common.cat021_df['mode_c_code'].iloc[index]
            alt_geo_ft = common.cat021_df['geometric_height'].iloc[index]

            assert not math.isnan(mc_ft)
            assert not math.isnan(alt_geo_ft)

            p = common.calc_p_from_baro_std(mc_ft * common.FT2M - self.geo_alt_correction_avg_m)
            geo_alt_rec_ft = common.calc_alt_from_estmations(p, self.est_p0, self.est_t0, self.est_alpha) * common.M2FT

            file1.write('{};{};{};{};{};{}\n'.format(
                common.time_str_from_seconds(common.cat021_df['time_of_day'].iloc[index]),
                mc_ft, mc_ft * common.FT2M, alt_geo_ft, geo_alt_rec_ft, p
            ))

        file1.close()

    def write(self, data_label_str, out_folder):

        assert self.finalized

        x1 = min(self.mcs_fts)
        x2 = max(self.mcs_fts)

        # calulate qnh curve
        baro_alts = np.arange(x1, x2 + 0.1, 0.1)

        if self.param_estimated:
            qnh_cor_alts_est = []

            for mc_ft in baro_alts:
                geo_alt_rec_ft = self.calc_alt_corr(mc_ft)
                qnh_cor_alts_est.append(geo_alt_rec_ft / 100)

            plt.plot([x / 100 for x in baro_alts], qnh_cor_alts_est, 'b-', label='P0 {:.2f} T0 {:.2f}'.format(
                self.est_p0, self.est_t0 - 273.15)) # alpha {:.5f} , self.est_alpha


        first_legend = plt.legend(loc="lower right")
        ax = plt.gca().add_artist(first_legend)

        scatter = plt.scatter(x=[x / 100 for x in self.mcs_fts], y=[x / 100 for x in self.alt_geos], c=self.calc_p0s, cmap='cool')

        plt.grid()
        plt.xlabel('Barometric Altitude [FL]')
        plt.ylabel('Geometric Altitude [FL]')
        fig = plt.gcf()
        fig.set_size_inches(common.plot_size_x, common.plot_size_y)

        common.create_path_if_required(out_folder)

        plt.savefig(out_folder + 'altitudes_cell_tod_{}_lat_{}_lon_{}.png'.format(
            common.time_str_from_seconds(self.tod_bins[self.tod_bin]),
            self.lat_bins[self.lat_bin], self.lon_bins[self.lon_bin]), dpi=common.plot_dpi, bbox_inches="tight")

        plt.close()

    def calc_alt_corr(self, mc_ft):

        assert self.finalized
        assert self.param_estimated

        # calculate pressure at mc for std formula wrt MSL (geoid)
        p = common.calc_p_from_baro_std(mc_ft * common.FT2M)

        # reconstruct geoidalt
        geo_alt_over_geoid = common.calc_alt_from_estmations(p, self.est_p0, self.est_t0, self.est_alpha)

        # remove goid alt to get ellipsoid alt
        return (geo_alt_over_geoid + self.geo_alt_correction_avg_m) * common.M2FT


