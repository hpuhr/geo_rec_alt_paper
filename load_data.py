#!/usr/bin/python3
import util.common as common
import pandas as pd

# common
# ds_id, uint
# time_of_day, float [s]
# aircraft_address, uint
# latitude, double, [deg]
# longitude, double, [deg]

#cat021
# geometric_height, float [ft]
# mode_c_code, float [ft]

#cat048
# mode_c_code, float [ft]
# mode_c_garbled, bool
# mode_c_valid, bool
# height_3d, float [ft], ref msl
# range, double [nm]
# azimuth, double [deg]

#cat062

def load_cat021(con):

    print("loading cat021 data")
    query_str = "SELECT time_of_day, aircraft_address, geometric_height, mode_c_code, latitude, longitude" \
                " from data_cat021 where mops_version = 2 AND nacp >= 4" \
                " AND mode_c_code IS NOT NULL AND geometric_height IS NOT NULL AND geometric_height != 32767*6.25" \
                " AND emitter_category in (3,5)" \
                " AND latitude >= {} AND latitude <= {}" \
                " AND longitude >= {} AND longitude <= {} ORDER BY time_of_day".format(
        common.latitude_min, common.latitude_max, common.longitude_min, common.longitude_max)

    print("query '{}'".format(query_str))
    df = pd.read_sql_query(query_str, con)

    # weird value Geometric Height 204793.75 = 32767*6.25

    print('got {} cat021 target reports'.format(df.shape[0]))
    return df

def load_cat048(con):

    print("loading cat048 data from {}".format(common.radar_ds_id))
    query_str = "SELECT time_of_day, aircraft_address, latitude, longitude, mode_c_code, height_3d, range, azimuth " \
                "from data_cat048 where ds_id = {} " \
                "AND mode_c_code IS NOT NULL AND mode_c_garbled = 0 AND mode_c_valid = 1 AND height_3d IS NOT NULL" \
                " AND aircraft_address IS NOT NULL AND mode_c_code > 10000 AND range < 62.5 ORDER BY time_of_day".format(common.radar_ds_id)

    print("query '{}'".format(query_str))
    df = pd.read_sql_query(query_str, con)

    print('got {} cat048 target reports'.format(df.shape[0]))
    print('latitude min/max {}/{}'.format(df['latitude'].min(), df['latitude'].max()))
    print('longitude min/max {}/{}'.format(df['longitude'].min(), df['longitude'].max()))

    return df

def load_cat062(con):

    print("loading cat062 data")
    query_str = "SELECT time_of_day, aircraft_address, latitude, longitude" \
                " from data_cat062 where aircraft_address IS NOT NULL AND latitude >= {} AND latitude <= {}" \
                " AND longitude >= {} AND longitude <= {} ORDER BY time_of_day".format(
        common.latitude_min, common.latitude_max, common.longitude_min, common.longitude_max)

    print("query '{}'".format(query_str))
    df = pd.read_sql_query(query_str, con)

    print('got {} cat062 target reports'.format(df.shape[0]))
    return df