import pandas as pd
from collections import OrderedDict
import gzip
import subprocess
import helper as h
import os
from time import time
import numpy as np
from collections import defaultdict
import sys
import argparse


def merge_maps(wave_maps, chunk):
    """
    Reads the SNPs from the .map files. Sorts them according to position keeping the original order in case of two variants with identical position.
    """

    def get_order(wave_maps):
        d = defaultdict(set)

        for index, wm in enumerate(wave_maps):
            grouped = wm.groupby('BP')
            for bp, group in grouped:
                if len(group) > 1:
                   for i, row in enumerate(group.iterrows()):
                       d[row[1].SNP] |= set(group.SNP[i:])

        d = {key: -len(val) for key, val in d.items()}
        return pd.DataFrame(d, index=['co']).transpose()

    order = get_order(wave_maps)

    merged_map = pd.concat(wave_maps)
    merged_map.drop_duplicates(inplace=True)
    merged_map = merged_map.reset_index(drop=True)

    order['SNP'] = order.index
    merged_map_order = pd.merge(order, merged_map, on='SNP', how='outer')

    merged_map_order.sort_values(by=['BP', 'co'], inplace=True)
    merged_map_order.reset_index(drop=True, inplace=True)
    merged_map_order.drop('co', axis=1, inplace=True)

    return merged_map_order


def write_dosage(d, s, include_inds, wave_fams):
    """
    Writes one row at a time to a merged dosage file, including only ind-columns in include_inds.
    """
    try:
        waves = []
        out = None
        out = open(s.mergepath, 'w')
        waves = [gzip.open(path + ".gz", 'rb') for path in s.wavepaths.filepaths]

        l_init = ["SNP", "A1", "A2"]
        l = []
        for wave in waves:
            l += wave.readline().strip().split()[3:]

        if include_inds:
            include_inds = set(include_inds)
            include_indices = []
            for i, elem in enumerate(l):
                if elem in include_inds:
                    include_indices.append(i - 1)
                    include_indices.append(i)
            l = [l[i] for i in include_indices]


        out.write(" ".join(l_init + l) + "\n")

        for snp, include_waves in d.items():
            l_init = [snp, "X", "Y"]
            l = []

            for index, wave in enumerate(waves):
                if index in include_waves:
                    line = wave.readline().strip().split()
                    l_init[1], l_init[2] = line[1], line[2]
                    for i in range(3, len(line), 2):
                        pair = line[i:i+2]
                        if not pair == ['NA','NA']:
                            p1, p2 = float(pair[0]), float(pair[1])
                            if p1 < genoprob and p2 < genoprob and 1- p1 - p2 < genoprob:
                                pair = ['NA', 'NA']
                        l += pair
                else:
                    l += [s.missing for _ in range(len(wave_fams[index])*2)]
            if include_inds:
                l = [l[i] for i in include_indices]

            out.write(" ".join(l_init + l) + "\n")

    finally:
        for wave in waves:
            if wave is not None:
                wave.close()

        if out is not None:
            out.close()


def make_dictionary(merged_map, s):
    """
    Creates a dictionary with snp ids as keys and the index of the waves that contain that snp as values.
    """
    d = OrderedDict((snp, set()) for snp in merged_map["SNP"])

    files = [gzip.open(path + ".gz",'rb') for path in s.wavepaths.filepaths]

    for i, f in enumerate(files):
        next(f)
        for line in f:
            line = line.split()
            d[line[0]].add(i)

    for f in files:
        f.close()

    return d


def merge_ngts(wave_ngts, merged_map):
    merged_ngts = merged_map
    ngtcols = ['GT' + str(i+1) for i in range(len(wave_ngts))]
    for i, ngt in enumerate(wave_ngts):
        ngt.rename(columns={'GT': 'GT' + str(i+1)}, inplace=True)
        merged_ngts = pd.merge(ngt, merged_ngts, on=['CHR', 'SNP', 'BP'], how='outer')
    merged_ngts = merged_ngts[['CHR', 'SNP','BP'] + ngtcols]
    merged_ngts['SUM'] = merged_ngts[ngtcols].sum(axis=1)
    merged_ngts['GT'] = np.where(merged_ngts[ngtcols].apply(lambda x: 0 in x.values, axis=1), 0, 1)
    return merged_ngts


def merge_infos(wave_infos, merged_map):
    merged_info = merged_map
    acols = ['a1_w' + str(i+1) for i in range(len(wave_ngts))]
    for i, info in enumerate(wave_infos):
        info.rename(columns={'a1': 'a1_w' + str(i+1), 'freq': 'freq_w' + str(i+1), 'pass': 'pass_w' + str(i+1),
                             'freq': 'freq_w' + str(i+1), 'info': 'info_w' + str(i+1)}, inplace=True)
        info.drop(['genotyped', 'a2'], inplace=True, axis=1)
        merged_info = pd.merge(info, merged_info, on=['CHR', 'SNP', 'BP'], how='outer')

    #assert all(len(set(row[1].tolist())) == 1 for row in merged_info[acols].iterrows())
    #assert merged_info[acols].apply( lambda x: len(x[-x.isnull()].unique()) == 1 , axis = 1).all() # Checks that all columns contains at most one variant (or NaN)
    return merged_info


s = h.parse()
s.wavepaths = h.replace(s)

genoprob = s.genoprob
s.missing = 'NA'
wave_maps = h.read_wave_maps(s.wavepaths)
wave_ngts = h.read_wave_ngts(s.wavepaths)
wave_infos = h.read_wave_infos(s.wavepaths)

merged_map = merge_maps(wave_maps, s.chunk)
merged_info = merge_infos(wave_infos, merged_map)

wave_fams = h.read_wave_fams(s.wavepaths)
merged_fam = pd.concat(wave_fams)

include_inds = h.read_include_inds(s.indlist)
if include_inds:
    merged_fam = merged_fam[merged_fam["indID"].isin(include_inds)]

d = make_dictionary(merged_map, s)

write_dosage(d, s, include_inds, wave_fams)

h.write_fam(s.mergepath, merged_fam)
h.write_map(s.mergepath, merged_map)


merged_ngts = merge_ngts(wave_ngts, merged_map)


h.write_info(s.mergepath, merged_info)
h.write_ngt(s.mergepath, merged_ngts)
h.write_ngt_matrix(s.mergepath, merged_ngts)

if os.path.exists(s.mergepath + ".gz"):
    os.remove(s.mergepath + ".gz")

subprocess.call(["gzip", s.mergepath])

