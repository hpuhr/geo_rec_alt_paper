import math
import os

DEG2RAD = 2*math.pi/360.0
RAD2DEG = 1.0/DEG2RAD
FT2M = 0.3048
M2FT = 1.0/FT2M

NM2M = 1852.0  # NM to meters
M2NM = 1.0/1852.0  # meters to NM
M_S2KNOTS = 3600.0 / 1852.0  # metres per second to knots

std_t0 = 273.15 + 15
std_p0 = 1013.25

plot_size_x = 6
plot_size_y = 4
plot_dpi = 300

geo_step = 0.25
tod_step = 60 * 60
num_bins = 40

cat021_df = None
cat048_df = None
cat062_df = None

latitude_bins = None
longitude_bins = None
time_of_day_bins = None

trk_chains = None # type: tracker_chains.TrackerChains
adsb_reconst_cells = None

## YOUR SETTINGS HERE

db_filename = "/yourpathtoyour/database.db"
egm96_filename = "/yourpathtoyour/egm96-5.pgm"
output_folder = 'output'

radar_ds_id = 1 # SAC * 255 + SIC

radar_lat = 4.1 # deg
radar_lon = 1.6 # deg
radar_alt = 42  # m over ellipsoid

# bouding rect deg
latitude_min = 1.0
latitude_max = 2.0
longitude_min = 1.0
longitude_max = 2.0

## DONE

def get_counts(data_column, bins=None):

    if bins is None:
        return data_column.value_counts().sort_index().to_frame(name="count").merge(
            data_column.value_counts(normalize=True).sort_index().to_frame(name="count(as %)"),
            left_index=True,
            right_index=True,
        )
    else:
        return data_column.value_counts(bins=bins).sort_index().to_frame(name="count").merge(
            data_column.value_counts(bins=bins, normalize=True).sort_index().to_frame(name="count(as %)"),
            left_index=True,
            right_index=True,
        )


def calc_p_from_baro_std(baro_m):
    return std_p0 * (1 - (0.0065*baro_m)/std_t0)**5.255 # calculate pressure at mc for std formula


def calc_p0_from_p_geo(p, geo_m):
    return p/((1 - (0.0065*geo_m)/std_t0)**5.255)


def calc_alt_from_p_p0(p, p0):
    return (std_t0/0.0065) * (1-(p/p0)**(1/5.255))

g0 = 9.81
R = 287

def calc_alt_from_estmations(p, p0, t0, alpha):
    return (t0/alpha) * (1-(p/p0)**((alpha * R) / g0))


def create_path_if_required(path_str):
    if not os.path.exists(path_str):
        os.makedirs(path_str)



def time_str_from_seconds(seconds, do_ms=True):
    hours = int(seconds / 3600)
    minutes = int((seconds-(hours*3600))/60)
    seconds_remain = seconds-(hours*3600)-minutes*60

    if do_ms:
        return "{0}:{1}:{2:06.3f}".format(str(int(hours)).zfill(2), str(int(minutes)).zfill(2), round(seconds_remain,3))
    else:
        return "{0}:{1}:{2}".format(str(int(hours)).zfill(2), str(int(minutes)).zfill(2),
                                    str(int(seconds_remain)).zfill(2))

def time_from_metar(value):
    hours = int(value / 100)
    minutes = int(value - (hours * 100))

    return hours * 3600 + minutes * 60
