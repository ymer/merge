import pandas as pd
from sys import argv
import argparse
import gzip
from itertools import chain
from os.path import basename
import helper as h
from time import time
import os
import numpy as np


def write_missing(s, wave_snps, merged_snps, wave_fams, merged_info):

    nwaves = len(wave_fams)

    pd.set_option('display.max_columns', None)
    l = []
    for snp in merged_snps:
        d = {}
        for i, wave in enumerate(s.wavepaths.path):
            wn = "info_w" + str(i+1)
            if snp in wave_snps[i]:
                if pd.isnull(merged_info[wn].loc[merged_info.SNP == snp]).any():
                    d[basename(wave)] = 0
                else:
                    d[basename(wave)] = 1
            else:
                d[basename(wave)] = 0
        d["SNP"] = snp
        l.append(d)
    df = pd.DataFrame(l)
    df.insert(1, "total", df.sum(axis=1))

    a = df["total"].value_counts(sort=True)
    a.name = "nsnps"
    a = a.reindex(range(nwaves,0,-1)).fillna(0).astype(int)
    a.to_csv(s.checkpath + "_missing_stats.txt", sep="\t", index_label=["present_in_waves"], header=True)

    df2 = df[df.total != nwaves]
    cols = ['SNP', 'total'] + [basename(wave) for wave in s.wavepaths.path]
    df2 = df2[cols]
    df2.to_csv(s.checkpath + "_missing.txt", sep="\t", index=False, na_rep="NA")

    for i, wave in enumerate(s.wavepaths.path):
        df[basename(wave)] = df[basename(wave)].apply(lambda x: x * len(wave_fams[i]))

    waves = range(1,nwaves + 1)
    for wave in waves:
        col = 'freq_w' + str(wave)
        merged_info[col] = np.where(merged_info[col] > 0.5, merged_info[col], 1 - merged_info[col])

    df3 = pd.merge(df, merged_info, on="SNP", how="inner")

    pd.set_option('display.max_columns', None)

    criteria = []
    for maf in [0, 0.01, 0.05]:
        for info in [0, 0.6, 0.8]:
            criteria.append(('INFO_' + str(info) + '_MAF_' + str(maf), info, maf))

    for c in criteria:
        colname, infomin, freqmin = c
        for wave in waves:
            wcol = basename(s.wavepaths.path[wave - 1])
            infocol = 'info_w' + str(wave)
            freqcol = 'freq_w' + str(wave)
            df3['w' + str(wave) + '_' + colname] = np.where((df3[infocol] > infomin) & (df3[freqcol] > freqmin), df3[wcol], 0)

        df3[colname] = df3[['w' + str(wave) + '_' + colname for wave in waves]].sum(axis=1)

    keep = ['SNP'] + [c[0] for c in criteria]
    df3 = df3[keep]

    df3.to_csv(s.checkpath + "_missing_with_indco.txt", sep="\t", index=False)


def original__inds_in_dosage_and_fam_are_identical():
    return all(inds == list(wave_fam["indID"]) for inds, wave_fam in zip(wave_inds, wave_fams))


def original__variants_in_dosage_and_map_are_identical():
    return all(snps == list(wave_map["SNP"])for snps, wave_map in zip(wave_snps, wave_maps))


def original__variants_are_sorted_by_position():
   return all(sorted(list(wave_map["BP"])) == list(wave_map["BP"]) for wave_map in wave_maps)


def original__indIDs_are_unique():
    comb = list(chain(*[list(wave_fam["indID"]) for wave_fam in wave_fams]))
    return len(comb) == len(set(comb))


def original__each_variant_has_just_one_position():
    df = wave_maps[0]
    df.rename(columns={"BP": "BP_0"}, inplace=True)
    for i, wave_map in enumerate(wave_maps[1:], start=1):
        wave_map.rename(columns={"BP": "BP_" + str(i)}, inplace=True)
        df = pd.merge(df, wave_map, how="outer", on="SNP")
    return df[["BP_" + str(i) for i in range(len(wave_maps))]].apply( lambda x: len(x[-x.isnull()].unique()) == 1 , axis = 1).all()


def merged__variants_are_sorted_by_position():
    return sorted(list(merged_map["BP"])) == list(merged_map["BP"])


def merged__inds_in_dosage_and_fam_are_identical():
    return merged_inds == list(merged_fam["indID"])


def merged__variants_in_dosage_and_map_are_identical():
    return merged_snps == list(merged_map["SNP"])


def compare__inds_in_merged_same_as_intersection_between_include_file_and_union_of_waves():
    return set(merged_inds) == set(chain(*wave_inds)) & set(include_inds)


def compare__inds_in_merged_same_as_in_union_of_waves():
    return set(merged_inds) == set(chain(*wave_inds))


def compare__variants_in_merged_same_as_union_of_variants_in_original():
    return merged_snps == list(merged_map["SNP"])


def compare__order_of_variants_in_waves_same_as_in_merged():

    for wave_index, wave_map in enumerate(wave_maps):
        wave_map["i"] = wave_map.index

        df = pd.merge(merged_map, wave_map, on="SNP", how="left").dropna()

        df['previ'] = df.i.shift()
        df2 = df[~(df.previ + 1 == df.i)]

        if not pd.Index(df.i).is_monotonic:
            return False

    return True


s = h.parse()
s.wavepaths = h.replace(s)

wave_maps = h.read_wave_maps(s.wavepaths)
wave_fams = h.read_wave_fams(s.wavepaths)
merged_map = h.read_map(s.mergepath)
merged_fam = h.read_fam(s.mergepath)

include_inds = h.read_include_inds(s.indlist)
if include_inds:
    merged_fam = merged_fam[merged_fam["indID"].isin(include_inds)]
wave_inds, wave_snps = h.read_wave_dosages(s.wavepaths.filepaths)
merged_inds, merged_snps = h.read_dosage(s.mergepath)
merged_info = h.read_info(s.mergepath)


checks = [original__inds_in_dosage_and_fam_are_identical,
          original__variants_in_dosage_and_map_are_identical,
          original__variants_are_sorted_by_position,
          original__indIDs_are_unique,
          original__each_variant_has_just_one_position,
          merged__variants_are_sorted_by_position,
          merged__inds_in_dosage_and_fam_are_identical,
          merged__variants_in_dosage_and_map_are_identical,
          compare__variants_in_merged_same_as_union_of_variants_in_original,
          compare__order_of_variants_in_waves_same_as_in_merged]

if s.indlist:
    checks += [compare__inds_in_merged_same_as_intersection_between_include_file_and_union_of_waves]
else:
    checks += [compare__inds_in_merged_same_as_in_union_of_waves]


with open(s.checkpath + "_check.txt", 'w') as out:
    for check in checks:
        result = check()
        out.write(str(result) + "\t" + check.__name__ + "\n")

write_missing(s, wave_snps, merged_snps, wave_fams, merged_info)

