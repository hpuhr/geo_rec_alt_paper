#!/usr/bin/python3
import util.common as common

from scipy.interpolate import interp1d
import bisect
import multiprocessing


def find_closest_value(alist, x):
    i = bisect.bisect_left(alist, x)
    if i >= len(alist):
        i = len(alist) - 1
    elif i and alist[i] - x > x - alist[i - 1]:
        i = i - 1
    return alist[i] # i

class TargetDataChain:
    def __init__(self, acad):
        self.acad = acad

        self.indexes = []

        self.tod_begin = None
        self.tod_end = None

        self.pos_tods = []
        self.positions = []

        self.f_pos_lin = None
        self.f_pos_quad = None
        self.f_pos_cub = None

        self.usable = False

    def add_index(self, index):

        self.indexes.append(index)

    def finalize(self):

        for index in self.indexes:

            assert common.cat062_df["aircraft_address"].iloc[index] == self.acad

            # same time issue
            if len(self.pos_tods) and common.cat062_df["time_of_day"].iloc[index] == self.pos_tods[-1]:
                continue

            self.pos_tods.append(common.cat062_df["time_of_day"].iloc[index])
            self.positions.append([common.cat062_df["latitude"].iloc[index],
                                   common.cat062_df["longitude"].iloc[index]])

        if len(self.pos_tods):
            self.tod_begin = self.pos_tods[0]
            self.tod_end = self.pos_tods[-1]

        if len(self.pos_tods) > 3:

            self.f_pos_lin = interp1d(self.pos_tods, self.positions, kind='linear', axis=0)

            self.f_pos_quad = interp1d(self.pos_tods, self.positions, kind='quadratic', axis=0)

            self.f_pos_cub = interp1d(self.pos_tods, self.positions, kind='cubic', axis=0)

            self.usable = True

        return self

    def exists_at(self, tod, max_delta):

        if not self.usable or not self.tod_begin or not self.tod_end:
            #print("nope2 usable {}".format(self.usable))
            return False

        if tod < self.tod_begin or tod > self.tod_end:
            #print("nope3")
            return False

        #print("nope4?")
        return abs(find_closest_value(self.pos_tods, tod) - tod) < max_delta

    def get_position_at(self, tod):

        #assert self.exists_at(tod)

        #assert self.f_pos_lin is not None
        #data_new2 = self.f_pos_lin(tod)

        #assert self.f_pos_quad is not None
        #data_new3 = self.f_pos_quad(tod)

        data_new4 = self.f_pos_cub(tod)

        lat = data_new4[0]
        lon = data_new4[1]

        # print('tod {} lat {} lon {}'.format(tod, lat, lon))

        return lat, lon

class TrackerChains:
    def __init__(self):

        self.chains = {}  # type: Dict[int, TargetDataChain]

    def load_data(self):

        no_acad_cnt = 0

        tmp_chains = {}

        for index in range(common.cat062_df.shape[0]):

            if index % 200000 == 0:
                print('creating tracker chains {}'.format(index))

            # data bins
            acad = common.cat062_df["aircraft_address"].iloc[index]

            if acad:
                if acad not in tmp_chains:
                    tmp_chains[acad] = TargetDataChain(acad)

                tmp_chains[acad].add_index(index)
            else:
                no_acad_cnt += 1

        print('creating {} radar chains done, {} of {} trs skipped'.format(
            len(tmp_chains), no_acad_cnt, common.cat062_df.shape[0]))

        print("finalizing chains")

        tmp_chain_vals = tmp_chains.values()

        pool = multiprocessing.Pool()  # initialise your pool
        finalized_chains = pool.map(TargetDataChain.finalize, tmp_chain_vals)
        pool.close()  # shut down the pool
        pool.join()

        num_usable = 0
        for chain in finalized_chains:
            self.chains[chain.acad] = chain

            if self.chains[chain.acad].usable:
                num_usable += 1

        print("finalizing {} chains done, usable {}".format(len(self.chains), num_usable))


    def has_position_at(self, acad, tod, time_delta):

        if acad not in self.chains:
            return False

        return self.chains[acad].exists_at(tod, time_delta)

    def get_position_at(self, acad, tod):

        assert acad in self.chains

        return self.chains[acad].get_position_at(tod)

def create():

    chains = TrackerChains()

    chains.load_data()

    return chains

