import pandas as pd
import argparse
from sys import argv
import gzip
import subprocess
import os

post = '.out.dosage'


def parse():
    parser = argparse.ArgumentParser()
    for arg in ['chunk', 'wavesfile', 'mergefolder', 'indlist', 'checkfolder', 'chunksfile']:
        parser.add_argument('--' + arg, dest=arg)

    for arg in ['genoprob']:
        parser.add_argument('--' + arg, dest=arg, type=float)

    s = parser.parse_args(argv[1:])

    if s.indlist == 'False':
        s.indlist = False

    if s.wavesfile and s.chunk:
        s.wavepaths = read_wavespaths(s.wavesfile, s.chunk)

    if s.mergefolder and s.chunk:
        s.mergepath = get_mergefolder(s)

    if s.checkfolder and s.chunk:
        s.checkpath = s.checkfolder + 'chunks/' + s.chunk

    return s


def get_mergefolder(s):
    return s.mergefolder + s.wavepaths.outname[0] + s.chunk + '.out.dosage'


def wavespaths_replace(wavespaths, replacefile):
    a = pd.read_csv(replacefile, sep='\s+')
    for i in range(len(a)):
        wavespaths.infopaths[a.wave[i]] = a.chunk[i]
        wavespaths.filepaths[a.wave[i]] = a.chunk[i]

    return wavespaths


def read_wavespaths(wavesfile, chunk):
    wavespaths = pd.read_csv(wavesfile, sep="\s+")
    wavespaths['filepaths'] = wavespaths.path.apply(lambda x: x + chunk + '.out.dosage')
    wavespaths['infopaths'] = wavespaths.infopath.apply(lambda x: x + chunk + '.out.dosage')
    wavespaths['outpaths'] = wavespaths.outname.apply(lambda x: x + chunk + '.out.dosage')
    return wavespaths

def replace(s):
    mf = s.mergefolder.split('/')[1]
    replacefile = 'temp/' + mf + '/' + s.chunk + '.txt'
    wavepaths = s.wavepaths
    if os.path.exists(replacefile):
        wavepaths = wavespaths_replace(s.wavepaths, replacefile)
    return wavepaths

def read_wave_fams(wavepaths, usecols=[]):
    return [read_fam(path, usecols) for path in wavepaths.filepaths]


def read_wave_maps(wavepaths):
    return [read_map(path) for path in wavepaths.filepaths]


def read_wave_ngts(wavepaths):
    return [read_ngt(path) for path in wavepaths.filepaths]


def read_ngt(path):
    return pd.read_csv(path + ".ngt", sep="\s+", header=None, names=["CHR", "SNP", "GEN", "BP", "GT"])


def read_wave_infos(wavepaths):
    return [read_info(path) for path in wavepaths.infopaths]


def read_info(path):
    info = pd.read_csv(path + ".info", sep="\s+")
    info = info.rename(columns={'POS': 'BP'})
    return info


def read_fam(path, usecols=[]):
    if usecols:
        return pd.read_csv(path + ".fam", sep="\s+", header=None, names= ["famID", "indID", "X", "Y", "gender", "aff"], usecols=usecols, dtype=str)
    else:
        return pd.read_csv(path + ".fam", sep="\s+", header=None, names= ["famID", "indID", "X", "Y", "gender", "aff"], dtype=str)


def read_map(path):
    return pd.read_csv(path + ".map", sep="\s+", header=None, usecols=[0,1,2,3], names=["CHR", "SNP", "GEN", "BP"])


def read_include_inds(include_inds_path):
    return list(pd.read_csv(include_inds_path, header=None, names=["famID", "indID", "pheno"], delimiter="\s+")["indID"]) if include_inds_path else None


def read_dosage(path):
    with gzip.open(path + ".gz",'rb') as inp:
        line = inp.readline().strip().split(" ")
        inds = line[4::2]
        snps = []
        for line in inp.readlines():
            snps.append(line.split()[0])

    return inds, snps


def read_wave_dosages(wavepaths):
    return zip(*[read_dosage(path) for path in wavepaths])


def write_fam(path, df):
    df.to_csv(path + ".fam", header=False, index=False, sep="\t")


def write_map(path, df):
    df = df[['CHR', 'SNP', 'GEN', 'BP']]
    df.to_csv(path + ".map", header=False, index=False, sep="\t")


def write_info(path, df):
    df.to_csv(path + '.info', sep=' ', index=False, na_rep='NA')


def read_assoc(path):
    return pd.read_csv(path + '.assoc.dosage', sep="\s+", dtype=object)


def read_chunks(path):
    return list(pd.read_csv(path, header=None)[0])


def extract_inds(infile, outfile, fam):

    inds = set(list(fam["famID"]) + list(fam["indID"]))

    with gzip.open(infile + '.gz', 'rb') as inp:
        line = inp.readline().split()

    cutlist = [0, 1, 2]
    for i, elem in enumerate(line):
        if elem in inds:
            cutlist.append(i)

    with open(outfile, 'w') as out:
        with gzip.open(infile + '.gz', 'rb') as inp:
            for line in inp:
                line = line.strip().split()
                useels = [line[i] for i in cutlist]
                if not set(useels[3:]) == set(["NA"]):
                    out.write(" ".join(useels) + "\n")

    subprocess.call(["gzip", outfile])


def write_ngt(path, df):
    df[['CHR', 'BP']] = df[['CHR', 'BP']].astype(int)
    df['GENO'] = 0
    df.to_csv(path + ".ngt", columns=['CHR', 'SNP', 'GENO', 'BP', 'GT'], header=False, index=False, sep=" ")


def write_ngt_matrix(path, df):
    df.drop(['CHR', 'BP', 'GENO'], axis=1, inplace=True)
    df.to_csv(path + ".ngt.matrix", index=False, sep=" ", na_rep='NaN', float_format='%1.0f')




