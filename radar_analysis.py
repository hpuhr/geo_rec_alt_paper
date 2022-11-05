#!/usr/bin/python3

import pymap3d
import numpy as np
import statistics
import matplotlib.pyplot as plt
import math
import pygeodesy
import multiprocessing
import itertools

import util.common as common

wgs84_ellispoid_semi_major = 6378137.0
wgs84_ellispoid_semi_minor = 6356752.314

def radius_at (latitude_deg):

    # see https://en.wikipedia.org/wiki/Reference_ellipsoid
    #  a and b are the equatorial radius (semi-major axis) and the polar radius (semi-minor axis),
    # N(lat) = a^2 / sqrt(a^2 * cos^2 (lat) + b^2 * sin^2 (lat))

    latitude_rad = latitude_deg * common.DEG2RAD

    return wgs84_ellispoid_semi_major ** 2 / math.sqrt(wgs84_ellispoid_semi_major ** 2 * math.cos(latitude_rad) ** 2
                                                       + wgs84_ellispoid_semi_minor ** 2 * math.sin(latitude_rad) ** 2)


def elev_angle(slant_range_m, altitude_m):
    return common.RAD2DEG * math.asin(
        ((altitude_m - common.radar_alt) - ((slant_range_m ** 2) / (2 * radius_at(common.radar_lat)))) / slant_range_m)

def lat_lon_from_polar_slant(azimuth_deg, s_range_m, altitude_m):

    elev_angle_deg = elev_angle(s_range_m, altitude_m)

    lat_calc, lon_calc, plot_h_calc = pymap3d.aer2geodetic(
        azimuth_deg, elev_angle_deg, s_range_m, common.radar_lat, common.radar_lon, common.radar_alt, deg=True)

    return lat_calc, lon_calc

def calculate_slice(index_list):

    calc_cnt = 0
    no_track_cnt = 0
    no_reconst_cnt = 0

    baro_cart_offset_slice = []
    baro_azm_offset_slice = []
    baro_range_offset_slice = []

    geom_cart_offset_slice = []
    geom_azm_offset_slice = []
    geom_range_offset_slice = []

    h3d_cart_offset_slice = []
    h3d_azm_offset_slice = []
    h3d_range_offset_slice = []

    for index in index_list:
        plt_tod = common.cat048_df["time_of_day"].iloc[index]
        plt_acad = common.cat048_df["aircraft_address"].iloc[index]
        plt_lat = common.cat048_df["latitude"].iloc[index]
        plt_lon = common.cat048_df["longitude"].iloc[index]
        plt_mc_m = common.cat048_df["mode_c_code"].iloc[index] * common.FT2M
        plot_slant_range_m = common.cat048_df["range"].iloc[index] * common.NM2M
        plot_azm_deg = (common.cat048_df["azimuth"].iloc[index])
        plot_height_3d_m = (common.cat048_df["height_3d"].iloc[index]) * common.FT2M

        if common.trk_chains.has_position_at(plt_acad, plt_tod, 6.0):
            trk_lat, trk_lon = common.trk_chains.get_position_at(plt_acad, plt_tod)

            # calc slant range corrected position with baro
            baro_elev_angle = elev_angle(plot_slant_range_m, plt_mc_m)

            baro_lat_calc, baro_lon_calc, baro_plot_h_calc = pymap3d.aer2geodetic(
                plot_azm_deg, baro_elev_angle, plot_slant_range_m, common.radar_lat, common.radar_lon,
                common.radar_alt, deg=True)

            # check if reconst cell available
            tod_bin = np.abs(np.array(common.time_of_day_bins) + common.tod_step / 2 - plt_tod).argmin()
            lat_bin = np.abs(np.array(common.latitude_bins) + common.geo_step / 2 - baro_lat_calc).argmin()
            lon_bin = np.abs(np.array(common.longitude_bins) + common.geo_step / 2 - baro_lon_calc).argmin()

            cell = common.adsb_reconst_cells[tod_bin][lat_bin][lon_bin]  # type: ReconstructorCell

            assert cell.finalized

            if not cell.param_estimated:
                no_reconst_cnt += 1
                continue

            plt_geoalt_m = (cell.calc_alt_corr(plt_mc_m * common.M2FT)) * common.FT2M

            if plt_geoalt_m >= plot_slant_range_m:
                print('impossible values sr {} geoalt {} baroalt {}'.format(
                    plot_slant_range_m, plt_geoalt_m, plt_mc_m))
                no_reconst_cnt += 1
                continue

            # calculate barometric offsets

            enu_calc_offset = pymap3d.geodetic2enu(baro_lat_calc, baro_lon_calc, 0,
                                                   trk_lat, trk_lon, 0, deg=True)
            calc_distance = math.sqrt(enu_calc_offset[0] ** 2 + enu_calc_offset[1] ** 2)

            baro_cart_offset_slice.append(calc_distance)

            trk_azimuth, trk_elevation_angle, trk_range = pymap3d.geodetic2aer(
                trk_lat, trk_lon, 0,
                common.radar_lat, common.radar_lon, 0, deg=True)

            angle_diff = trk_azimuth - plot_azm_deg
            angle_diff = (angle_diff + 180) % 360 - 180

            baro_azm_offset_slice.append(angle_diff)
            geom_azm_offset_slice.append(angle_diff)

            tmp, tmp, rad_range_cor_baro = pymap3d.geodetic2aer(
                baro_lat_calc, baro_lon_calc, 0, common.radar_lat, common.radar_lon, 0, deg=True)
            baro_range_offset_slice.append(trk_range - rad_range_cor_baro)

            # calculate geometric offsets

            geo_elev_angle = elev_angle(plot_slant_range_m, plt_geoalt_m)

            geo_lat_calc, geo_lon_calc, geo_plot_h_calc = pymap3d.aer2geodetic(
                plot_azm_deg, geo_elev_angle, plot_slant_range_m, common.radar_lat, common.radar_lon,
                common.radar_alt, deg=True)

            enu_calc_offset = pymap3d.geodetic2enu(geo_lat_calc, geo_lon_calc, 0,
                                                   trk_lat, trk_lon, 0, deg=True)
            calc_distance = math.sqrt(enu_calc_offset[0] ** 2 + enu_calc_offset[1] ** 2)

            geom_cart_offset_slice.append(calc_distance)

            rad_azimuth_cor, rad_elevation_angle_cor, rad_range_cor = pymap3d.geodetic2aer(
                geo_lat_calc, geo_lon_calc, 0,
                common.radar_lat, common.radar_lon, 0, deg=True)

            angle_diff = trk_azimuth - rad_azimuth_cor
            angle_diff = (angle_diff + 180) % 360 - 180

            geom_azm_offset_slice.append(angle_diff)
            geom_range_offset_slice.append(trk_range - rad_range_cor)

            # calculate height 3d offsets

            h3d_elev_angle = elev_angle(
                plot_slant_range_m, plot_height_3d_m - common.cat021_df['geometric_height_correction_m'].iloc[index])
            # correct for msl

            h3d_lat_calc, h3d_lon_calc, h3d_plot_h_calc = pymap3d.aer2geodetic(
                plot_azm_deg, h3d_elev_angle, plot_slant_range_m, common.radar_lat, common.radar_lon,
                common.radar_alt, deg=True)

            enu_calc_offset = pymap3d.geodetic2enu(h3d_lat_calc, h3d_lon_calc, 0,
                                                   trk_lat, trk_lon, 0, deg=True)
            calc_distance = math.sqrt(enu_calc_offset[0] ** 2 + enu_calc_offset[1] ** 2)

            h3d_cart_offset_slice.append(calc_distance)

            rad_azimuth_cor, rad_elevation_angle_cor, rad_range_cor = pymap3d.geodetic2aer(
                h3d_lat_calc, h3d_lon_calc, 0,
                common.radar_lat, common.radar_lon, 0, deg=True)

            angle_diff = trk_azimuth - rad_azimuth_cor
            angle_diff = (angle_diff + 180) % 360 - 180

            h3d_azm_offset_slice.append(angle_diff)
            h3d_range_offset_slice.append(trk_range - rad_range_cor)

            calc_cnt += 1

            # print("polar offset azm {} range {}".format(trk_azimuth-plt_azimuth, trk_range-plt_range))
        else:
            no_track_cnt += 1

    return calc_cnt, no_track_cnt, no_reconst_cnt, \
           baro_cart_offset_slice, baro_azm_offset_slice, baro_range_offset_slice, \
           geom_cart_offset_slice, geom_azm_offset_slice, geom_range_offset_slice, \
           h3d_cart_offset_slice, h3d_azm_offset_slice, h3d_range_offset_slice

def calc_original_offsets():

    baro_cart_offsets = {}
    baro_azm_offsets = {}
    baro_range_offsets = {}

    geom_cart_offsets = {}
    geom_azm_offsets = {}
    geom_range_offsets = {}

    h3d_cart_offsets = {}
    h3d_azm_offsets = {}
    h3d_range_offsets = {}

    no_track_cnt = 0
    no_reconst_cnt = 0
    calc_cnt = 0

    tod_slices = {} # tod slice -> index list

    print("calculating radar offsets")

    # calculate time slices

    for index in range(common.cat048_df.shape[0]):
        plt_tod = common.cat048_df["time_of_day"].iloc[index]

        tod_slice = int(plt_tod / (30 * 60)) * 30 * 60

        if tod_slice not in tod_slices:
            tod_slices[tod_slice] = []

        tod_slices[tod_slice].append(index)

    tmp_slices = tod_slices.values()
    pool = multiprocessing.Pool()  # initialise your pool
    results = pool.map(calculate_slice, tmp_slices)
    pool.close()  # shut down the pool
    pool.join()

    assert len(results) == len(tod_slices)

    result_cnt = 0

    for tod_slice, indexes in tod_slices.items():
        calc_cnt_slice, no_track_cnt_slice, no_reconst_cnt_slice, \
        baro_cart_offset_slice, baro_azm_offset_slice, baro_range_offset_slice, \
        geom_cart_offset_slice, geom_azm_offset_slice, geom_range_offset_slice, \
           h3d_cart_offset_slice, h3d_azm_offset_slice, h3d_range_offset_slice = results[result_cnt]

        no_track_cnt += no_track_cnt_slice
        no_reconst_cnt += no_reconst_cnt_slice
        calc_cnt += calc_cnt_slice

        baro_cart_offsets[tod_slice] = baro_cart_offset_slice
        baro_azm_offsets[tod_slice] = baro_azm_offset_slice
        baro_range_offsets[tod_slice] = baro_range_offset_slice

        geom_cart_offsets[tod_slice] = geom_cart_offset_slice
        geom_azm_offsets[tod_slice] = geom_azm_offset_slice
        geom_range_offsets[tod_slice] = geom_range_offset_slice

        h3d_cart_offsets[tod_slice] = h3d_cart_offset_slice
        h3d_azm_offsets[tod_slice] = h3d_azm_offset_slice
        h3d_range_offsets[tod_slice] = h3d_range_offset_slice

        result_cnt += 1


    print('calculated {} radar offsets of {}, track misses {} reconst misses {}'.format(
        calc_cnt, common.cat048_df.shape[0], no_track_cnt, no_reconst_cnt))

    # print cart

    tods = []
    baro_values = []
    baro_stddevs = []
    geom_values = []
    geom_stddevs = []
    h3d_values = []
    h3d_stddevs = []

    for tod_slice, offset_list in baro_cart_offsets.items():

        if len(offset_list) < 2:
            continue

        tods.append(tod_slice)
        baro_values.append(statistics.mean(offset_list))
        baro_stddevs.append(statistics.stdev(offset_list))

    for tod_slice, offset_list in geom_cart_offsets.items():

        if len(offset_list) < 2:
            continue

        geom_values.append(statistics.mean(offset_list))
        geom_stddevs.append(statistics.stdev(offset_list))

    for tod_slice, offset_list in h3d_cart_offsets.items():

        if len(offset_list) < 2:
            continue

        h3d_values.append(statistics.mean(offset_list))
        h3d_stddevs.append(statistics.stdev(offset_list))

    baro_values_full = list(itertools.chain.from_iterable(baro_cart_offsets.values()))
    geom_values_full = list(itertools.chain.from_iterable(geom_cart_offsets.values()))
    h3d_values_full = list(itertools.chain.from_iterable(h3d_cart_offsets.values()))


    plt.plot(tods, h3d_values, 'go-', label='H.3D Avg. {:.2f}'.format(statistics.mean(h3d_values_full)))
    plt.plot(tods, baro_values, 'ro-', label='Baro. Avg. {:.2f}'.format(statistics.mean(baro_values_full)))
    plt.plot(tods, geom_values, 'bo-', label='Geom. Avg. {:.2f}'.format(statistics.mean(geom_values_full)))
    plt.errorbar(tods, baro_values, yerr=baro_stddevs, fmt='r.', capsize=2, elinewidth=0, alpha=0.5)
    plt.errorbar(tods, geom_values, yerr=geom_stddevs, fmt='b.', capsize=2, elinewidth=0, alpha=0.5)

    plt.legend(loc='best')
    plt.xlabel('Time [s]')
    plt.ylabel('Cartesian Offset [m]')

    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.gca().set_xticklabels([common.time_str_from_seconds(x, False) for x in plt.gca().get_xticks()])
    plt.savefig(common.output_folder + '/' + 'radar_cart_offset', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()

    print("cart offset baro avg {} stddev {}".format(statistics.mean(
        baro_values_full), statistics.stdev(baro_values_full)))
    print("cart offset geom avg {} stddev {}".format(
        statistics.mean(geom_values_full), statistics.stdev(geom_values_full)))
    print("cart offset h3d avg {} stddev {}".format(
        statistics.mean(h3d_values_full), statistics.stdev(h3d_values_full)))

    # print azm

    tods = []
    baro_values = []
    baro_stddevs = []
    geom_values = []
    geom_stddevs = []
    h3d_values = []
    h3d_stddevs = []

    for tod_slice, offset_list in baro_azm_offsets.items():

        if len(offset_list) < 2:
            continue

        #print('baro azm offset tod {} avg {} std.dev. {} len {}'.format(
        #    tod_slice, statistics.mean(offset_list), statistics.stdev(offset_list), len(offset_list)))
        tods.append(tod_slice)
        baro_values.append(statistics.mean(offset_list))
        baro_stddevs.append(statistics.stdev(offset_list))

    for tod_slice, offset_list in geom_azm_offsets.items():

        if len(offset_list) < 2:
            continue

        #print('geom azm offset tod {} avg {} std.dev. {} len {}'.format(
        #    tod_slice, statistics.mean(offset_list), statistics.stdev(offset_list), len(offset_list)))
        geom_values.append(statistics.mean(offset_list))
        geom_stddevs.append(statistics.stdev(offset_list))

    for tod_slice, offset_list in h3d_azm_offsets.items():

        if len(offset_list) < 2:
            continue

        h3d_values.append(statistics.mean(offset_list))
        h3d_stddevs.append(statistics.stdev(offset_list))

    baro_values_full = list(itertools.chain.from_iterable(baro_azm_offsets.values()))
    geom_values_full = list(itertools.chain.from_iterable(geom_azm_offsets.values()))
    h3d_values_full = list(itertools.chain.from_iterable(h3d_azm_offsets.values()))

    plt.plot(tods, h3d_values, 'go-', label='H.3D Avg. {:.2f}'.format(statistics.mean(h3d_values_full)))
    plt.plot(tods, baro_values, 'ro-', label='Baro. Avg. {:.2f}'.format(statistics.mean(baro_values_full)))
    plt.plot(tods, geom_values, 'bo-', label='Geom. Avg. {:.2f}'.format(statistics.mean(geom_values_full)))
    plt.errorbar(tods, baro_values, yerr=baro_stddevs, fmt='r.', capsize=2, elinewidth=0, alpha=0.5)
    plt.errorbar(tods, geom_values, yerr=geom_stddevs, fmt='b.', capsize=2, elinewidth=0, alpha=0.5)

    plt.legend(loc='best')
    plt.xlabel('Time [s]')
    plt.ylabel('Azimuth Offset [deg]')

    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.gca().set_xticklabels([common.time_str_from_seconds(x, False) for x in plt.gca().get_xticks()])
    plt.savefig(common.output_folder + '/' + 'radar_azm_offset', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()

    print("azm offset baro avg {} stddev {}".format(statistics.mean(
        baro_values_full), statistics.stdev(baro_values_full)))
    print("azm offset geom avg {} stddev {}".format(
        statistics.mean(geom_values_full), statistics.stdev(geom_values_full)))
    print("azm offset h3d avg {} stddev {}".format(
        statistics.mean(h3d_values_full), statistics.stdev(h3d_values_full)))

    # print range

    tods = []
    baro_values = []
    baro_stddevs = []
    geom_values = []
    geom_stddevs = []
    h3d_values = []
    h3d_stddevs = []

    for tod_slice, offset_list in baro_range_offsets.items():

        if len(offset_list) < 2:
            continue

        #print('baro range offset tod {} avg {} std.dev. {} len {}'.format(
        #    tod_slice, statistics.mean(offset_list), statistics.stdev(offset_list), len(offset_list)))
        tods.append(tod_slice)
        baro_values.append(statistics.mean(offset_list))
        baro_stddevs.append(statistics.stdev(offset_list))

    for tod_slice, offset_list in geom_range_offsets.items():

        if len(offset_list) < 2:
            continue

        #print('geom range offset tod {} avg {} std.dev. {} len {}'.format(
        #    tod_slice, statistics.mean(offset_list), statistics.stdev(offset_list), len(offset_list)))
        geom_values.append(statistics.mean(offset_list))
        geom_stddevs.append(statistics.stdev(offset_list))

    for tod_slice, offset_list in h3d_range_offsets.items():

        if len(offset_list) < 2:
            continue

        h3d_values.append(statistics.mean(offset_list))
        h3d_stddevs.append(statistics.stdev(offset_list))

    baro_values_full = list(itertools.chain.from_iterable(baro_range_offsets.values()))
    geom_values_full = list(itertools.chain.from_iterable(geom_range_offsets.values()))
    h3d_values_full = list(itertools.chain.from_iterable(h3d_range_offsets.values()))

    plt.plot(tods, h3d_values, 'go-', label='H.3D Avg. {:.2f}'.format(statistics.mean(h3d_values_full)))
    plt.plot(tods, baro_values, 'ro-', label='Baro. Avg. {:.2f}'.format(statistics.mean(baro_values_full)))
    plt.plot(tods, geom_values, 'bo-', label='Geom. Avg. {:.2f}'.format(statistics.mean(geom_values_full)))
    plt.errorbar(tods, baro_values, yerr=baro_stddevs, fmt='r.', capsize=2, elinewidth=0, alpha=0.5)
    plt.errorbar(tods, geom_values, yerr=geom_stddevs, fmt='b.', capsize=2, elinewidth=0, alpha=0.5)

    plt.legend(loc='best')
    plt.xlabel('Time [s]')
    plt.ylabel('Range Offset [m]')
    fig = plt.gcf()
    fig.set_size_inches(common.plot_size_x, common.plot_size_y)
    plt.gca().set_xticklabels([common.time_str_from_seconds(x, False) for x in plt.gca().get_xticks()])
    plt.savefig(common.output_folder + '/' + 'radar_range_offset', dpi=common.plot_dpi, bbox_inches="tight")
    plt.close()

    print("range offset baro avg {} stddev {}".format(statistics.mean(
        baro_values_full), statistics.stdev(baro_values_full)))
    print("range offset geom avg {} stddev {}".format(
        statistics.mean(geom_values_full), statistics.stdev(geom_values_full)))
    print("range offset h3d avg {} stddev {}".format(
        statistics.mean(h3d_values_full), statistics.stdev(h3d_values_full)))


def test_coordinates():

    #trk_diffs = []
    azm_diffs = []
    rng_diffs = []
    elv_diffs = []
    cart_offs = []
    cart_offs_e = []
    cart_offs_n = []

    for index in range(common.cat048_df.shape[0]):
        plt_tod = common.cat048_df["time_of_day"].iloc[index]
        plt_acad = common.cat048_df["aircraft_address"].iloc[index]
        plt_lat = common.cat048_df["latitude"].iloc[index]
        plt_lon = common.cat048_df["longitude"].iloc[index]
        plt_mc_m = common.cat048_df["mode_c_code"].iloc[index] * common.FT2M
        plot_slant_range_m = common.cat048_df["range"].iloc[index] * common.NM2M
        plot_azm_deg = (common.cat048_df["azimuth"].iloc[index])
        plot_height_3d_m = (common.cat048_df["height_3d"].iloc[index]) * common.FT2M

        if not common.trk_chains.has_position_at(plt_acad, plt_tod, 6.0):
            continue

        trk_lat, trk_lon = common.trk_chains.get_position_at(plt_acad, plt_tod)

        # check if reconst cell available
        tod_bin = np.abs(np.array(common.time_of_day_bins) + common.tod_step / 2 - plt_tod).argmin()
        lat_bin = np.abs(np.array(common.latitude_bins) + common.geo_step / 2 - plt_lat).argmin()
        lon_bin = np.abs(np.array(common.longitude_bins) + common.geo_step / 2 - plt_lon).argmin()

        cell = common.adsb_reconst_cells[tod_bin][lat_bin][lon_bin]  # type: ReconstructorCell

        assert cell.finalized

        if not cell.param_estimated:
            continue

        plt_geoalt_m = (cell.calc_alt_corr(plt_mc_m * common.M2FT)) * common.FT2M

        if plt_geoalt_m >= plot_slant_range_m:
            continue

        #print("baro {} geo {} diff {}".format(plt_mc_m, plt_geoalt_m, plt_mc_m-plt_geoalt_m))

        geo_elev_angle = elev_angle(plot_slant_range_m, plt_geoalt_m)

        lat_calc, lon_calc, plot_h_calc = pymap3d.aer2geodetic(
            plot_azm_deg, geo_elev_angle, plot_slant_range_m, common.radar_lat, common.radar_lon, common.radar_alt, deg=True)

        plot_azimuth_calc, plot_elevation_angle_calc, plot_range_calc = pymap3d.geodetic2aer(
            lat_calc, lon_calc, 0,
            common.radar_lat, common.radar_lon, 0, deg=True)

        trk_azimuth_calc, trk_elevation_angle_calc, trk_range_calc = pymap3d.geodetic2aer(
            trk_lat, trk_lon, 0,
            common.radar_lat, common.radar_lon, 0, deg=True)

        elv_diffs.append(plot_elevation_angle_calc - trk_elevation_angle_calc)
        #print('azm rad {} calc {} diff {}'.format(plot_azm_deg, rad_azimuth_calc, plot_azm_deg - rad_azimuth_calc))

        angle_diff = plot_azm_deg - trk_azimuth_calc
        angle_diff = (angle_diff + 180) % 360 - 180

        azm_diffs.append(angle_diff)

        #print('rng rad {} calc {} diff {}'.format(plot_slant_range_m, rad_range_calc,
        #                                          plot_slant_range_m - rad_range_calc))
        rng_diffs.append(plot_range_calc - trk_range_calc)
        #print('elev geo {} calc {} diff {}'.format(geo_elev_angle, rad_elevation_angle_calc,
        #                                          geo_elev_angle - rad_elevation_angle_calc))

        enu_calc_offset = pymap3d.geodetic2enu(trk_lat, trk_lon, 0,
                                               lat_calc, lon_calc, 0, deg=True)
        calc_distance = math.sqrt(enu_calc_offset[0] ** 2 + enu_calc_offset[1] ** 2)

        cart_offs.append(calc_distance)

        cart_offs_e.append(enu_calc_offset[0])
        cart_offs_n.append(enu_calc_offset[1])

    print('azm diff avg {} stddev {}'.format(statistics.mean(azm_diffs), statistics.stdev(azm_diffs)))
    print('rng diff avg {} stddev {}'.format(statistics.mean(rng_diffs), statistics.stdev(rng_diffs)))
    print('elv diff avg {} stddev {}'.format(statistics.mean(elv_diffs), statistics.stdev(elv_diffs)))
    print('cart_offs avg {} stddev {}'.format(statistics.mean(cart_offs), statistics.stdev(cart_offs)))
    print('cart_offs_e avg {} stddev {}'.format(statistics.mean(cart_offs_e), statistics.stdev(cart_offs_e)))
    print('cart_offs_n avg {} stddev {}'.format(statistics.mean(cart_offs_n), statistics.stdev(cart_offs_n)))

    lat, lon, alt = pymap3d.enu2geodetic(statistics.mean(cart_offs_e), statistics.mean(cart_offs_n), 0,
                                         common.radar_lat, common.radar_lon, 0, deg=True)
    print('new coords lat {} lon {}'.format(lat, lon))

