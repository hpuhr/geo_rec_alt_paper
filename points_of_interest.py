#!/usr/bin/python3

import pygeodesy
import numpy as np
import matplotlib.pyplot as plt

import util.common as common

def calculate():

    ginterpolator = pygeodesy.GeoidPGM(common.egm96_filename)

    points_of_interest = {'loww': (48.1126, 16.5755)}

    # time, t0 air, t0 dew, qnh
    ref_metar = {'loww': [[common.time_from_metar(20), 21, 17, 1011],
                          [common.time_from_metar(50), 21, 17, 1011],
                          [common.time_from_metar(120), 21, 17, 1011],
                          [common.time_from_metar(150), 20, 17, 1011],
                          [common.time_from_metar(220), 20, 16, 1011],
                          [common.time_from_metar(250), 20, 16, 1011],
                          [common.time_from_metar(320), 19, 16, 1011],
                          [common.time_from_metar(350), 19, 16, 1011],
                          [common.time_from_metar(420), 19, 16, 1012],
                          [common.time_from_metar(450), 19, 16, 1012],
                          [common.time_from_metar(520), 19, 16, 1012],
                          [common.time_from_metar(550), 20, 17, 1012],
                          [common.time_from_metar(620), 20, 17, 1012],
                          [common.time_from_metar(650), 22, 18, 1012],
                          [common.time_from_metar(720), 24, 18, 1013],
                          [common.time_from_metar(750), 24, 17, 1013],
                          [common.time_from_metar(820), 25, 17, 1013],
                          [common.time_from_metar(850), 25, 17, 1013],
                          [common.time_from_metar(920), 26, 17, 1013],
                          [common.time_from_metar(950), 27, 16, 1012],
                          [common.time_from_metar(1020), 27, 16, 1012],
                          [common.time_from_metar(1050), 27, 16, 1012],
                          [common.time_from_metar(1120), 28, 15, 1012],
                          [common.time_from_metar(1150), 28, 16, 1012],
                          [common.time_from_metar(1220), 29, 17, 1012],
                          [common.time_from_metar(1250), 28, 16, 1012],
                          [common.time_from_metar(1320), 28, 16, 1011],
                          [common.time_from_metar(1350), 28, 16, 1011],
                          [common.time_from_metar(1420), 28, 16, 1011],
                          [common.time_from_metar(1450), 28, 16, 1011],
                          [common.time_from_metar(1520), 28, 17, 1011],
                          [common.time_from_metar(1550), 29, 17, 1011],
                          [common.time_from_metar(1620), 27, 17, 1011],
                          [common.time_from_metar(1650), 27, 17, 1011],
                          [common.time_from_metar(1720), 26, 17, 1012],
                          [common.time_from_metar(1750), 25, 17, 1012],
                          [common.time_from_metar(1820), 23, 18, 1012],
                          [common.time_from_metar(1850), 23, 18, 1012],
                          [common.time_from_metar(1920), 23, 18, 1012],
                          [common.time_from_metar(1950), 24, 18, 1012],
                          [common.time_from_metar(2020), 23, 17, 1013],
                          [common.time_from_metar(2050), 22, 16, 1013],
                          [common.time_from_metar(2120), 21, 15, 1013],
                          [common.time_from_metar(2150), 20, 16, 1013],
                          [common.time_from_metar(2220), 20, 16, 1013],
                          [common.time_from_metar(2250), 20, 16, 1014],
                          [common.time_from_metar(2320), 19, 16, 1014],
                          [common.time_from_metar(2350), 18, 17, 1014],
                          [common.time_from_metar(2350), 18, 17, 1014]
                          ]}

    for p_lbl, latlong in points_of_interest.items():

        ref_times = [i[0] for i in ref_metar[p_lbl]]
        ref_t0s = [i[1] for i in ref_metar[p_lbl]]
        ref_p0s = [i[3] for i in ref_metar[p_lbl]]

        times = []
        p0s = []
        t0s = []

        lat, lon = latlong

        geoid_h = ginterpolator.height(lat, lon)

        print('doing {} lat {} lon {} h {}'.format(p_lbl, lat, lon, geoid_h))

        #self.lat_min, self.lat_max, self.lat_bins, self.lat_data_bins

        lat_bin = np.abs(np.array(common.latitude_bins) + common.geo_step/2 - lat).argmin()
        lon_bin = np.abs(np.array(common.longitude_bins) + common.geo_step/2 - lon).argmin()

        print('lat bin {} val {} lon bin {} val {}'.format(
            lat_bin, common.latitude_bins[lat_bin], lon_bin, common.longitude_bins[lon_bin]))

        for time_cells in common.adsb_reconst_cells:
            cell = time_cells[lat_bin][lon_bin] # type: ReconstructorCell

            assert cell.finalized

            if cell.param_estimated:
                #cell.do_param_estimation()

                #print('lat {} {} lon {} {}'.format(cell.lat_min, cell.lat_max, cell.lon_min, cell.lon_max))

                print('time {} count {} p0 {} T0 {} alpha {}'.format(
                    common.time_str_from_seconds(cell.tod_bins[cell.tod_bin]),
                    len(cell.indexes), cell.est_p0,
                    cell.est_t0 - 273.15, cell.est_alpha)) # cell.do_param_estimation()

                times.append(cell.tod_bins[cell.tod_bin] + common.tod_step / 2)
                p0s.append(cell.est_p0)
                t0s.append(cell.est_t0 - 273.15)

                #cell.write_calc('reconst', common.output_folder + '/' + p_lbl + '/')
                cell.write('reconst', common.output_folder + '/' + p_lbl + '/')
            else:
                print("time {} count {} no result".format(common.time_str_from_seconds(cell.tod_bins[cell.tod_bin]),
                                                          len(cell.indexes)))

        if len(times):

            #print("times {}: {}".format(len(times), times))
            #print("p0s {}: {}".format(len(p0s), p0s))
            #print("t0s {}: {}".format(len(t0s), t0s))

            #write p0s
            plt.plot(times, p0s, 'bo-', label='Reconst')
            plt.plot(ref_times, ref_p0s, 'ko-', label='MET')

            plt.grid()
            plt.xlabel('Time')
            plt.ylabel('P0')
            plt.legend(loc="lower right")
            plt.gca().set_xticklabels([common.time_str_from_seconds(x, False) for x in plt.gca().get_xticks()])
            fig = plt.gcf()
            fig.set_size_inches(common.plot_size_x, common.plot_size_y)

            plt.savefig(common.output_folder + '/loww_p0.png', dpi=common.plot_dpi, bbox_inches="tight")
            plt.close()

            # write t0s
            plt.plot(times, t0s, 'bo-', label='Reconst')
            plt.plot(ref_times, ref_t0s, 'ko-', label='MET')

            plt.grid()
            plt.xlabel('Time')
            plt.ylabel('T0')
            plt.legend(loc="lower right")
            plt.gca().set_xticklabels([common.time_str_from_seconds(x, False) for x in plt.gca().get_xticks()])
            fig = plt.gcf()
            fig.set_size_inches(common.plot_size_x, common.plot_size_y)

            plt.savefig(common.output_folder + '/loww_t0.png', dpi=common.plot_dpi, bbox_inches="tight")
            plt.close()

