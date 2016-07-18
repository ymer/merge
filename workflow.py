from gwf import *
import pandas as pd
from os import path, mkdir
import random
import itertools as it
import time
from shutil import copyfile
import pickle


def prepare():
    o = {'output': 'prepare'}
    cmd = "python scripts/prepare.py --chunksfile {} --wavesfile {} --mergefolder {}".format(chunksfile, wavesfile, mergefolder)
    return o, cmd


def merge(chunk):
    o = {'input': 'prepare', 'output': 'merge_' + chunk}
    cmd = "python scripts/merge.py --chunk {} --indlist {} --wavesfile {} --mergefolder {} --genoprob {}".format(chunk, s.indlist, wavesfile, mergefolder, s.genoprob)
    return o, cmd


def check(chunk):
    o = {'input': 'merge_' + chunk, 'output': 'check_' + chunk}
    cmd = "python scripts/check.py --chunk {} --indlist {} --wavesfile {} --mergefolder {} --checkfolder {}".format(chunk, s.indlist, wavesfile, mergefolder, checkfolder)
    return o, cmd


def collect_check():

    o = {'input': ['check_' + chunk for chunk in chunks], 'output': 'collect_check'}
    cmd = "python scripts/collect_check.py --chunksfile {} --checkfolder {}".format(chunksfile, checkfolder)
    return o, cmd


def write_covar():
    o = {'input': ['merge_' + chunk for chunk in chunks], 'output': 'covar'}
    cmd = "python scripts/write_covar.py --chunk {} --checkfolder {} --wavesfile {} --mergefolder {}".format(chunk, checkfolder, wavesfile, mergefolder)
    return o, cmd


def create_dirs():
    def make_path(pth):
        if not path.exists(pth):
            mkdir(pth)

    for pth in ['output', 'output/' + folder, mergefolder, checkfolder, checkfolder + 'chunks']:
        make_path(pth)


s = pickle.load(open('files/config.pl', 'rb'))

folder = s.vers + '_' + s.date + '_all_' + str(s.genoprob) + '/'
chunksfile = 'files/chunks' + s.chunks + '.txt'
checkfolder = 'output/' + folder + 'check/'
mergefolder = 'output/' + folder + 'merged/'
wavesfile = 'files/file_locations_' + s.vers + '.txt'
create_dirs()

chunks = list(pd.read_csv(chunksfile, header=None)[0])

w = 0
if w:
    o, cmd = merge(chunks[0])
    print o
    print
    print cmd
else:
    copyfile('files/readme.txt', 'output/' + folder + 'readme.txt')
    target('Prepare') << prepare()
    for chunk in chunks:
        target('Merge_' + chunk, walltime='10:00:00') << merge(chunk)
        if s.run_check:
            target('Check_' + chunk, memory="30g", walltime='10:00:00') << check(chunk)

    if s.run_check:
        target('Collect_check', memory='20g') << collect_check()
    target('Write_covar') << write_covar()

