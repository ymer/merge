import helper as h
from os import path
import pandas as pd
import numpy as np
import os

s = h.parse()

chunks = list(pd.read_csv(s.chunksfile, header=None)[0])
wp0 = h.read_wavespaths(s.wavesfile, chunks[0])
mf = s.mergefolder.split('/')[1]
if not os.path.exists('temp/' + mf):
    os.makedirs('temp/' + mf)



def make_nan(wavespaths, chunk, i):

    chunkfam = h.read_fam(wp0.filepaths[i])
    chunkmap = h.read_map(wavespaths.filepaths[0])
    chunkngt = h.read_ngt(wavespaths.filepaths[0])
    chunkinfo = h.read_info(wavespaths.infopaths[0])
    chunkinfo['info'] = "NA"
    chunkinfo['freq'] = "NA"
    chunkinfo['a1'] = "NA"
    chunkinfo['a2'] = "NA"

    dosage = pd.read_csv(wavespaths.filepaths[0] + '.gz', sep=' ', compression='gzip', usecols=['SNP', 'A1', 'A2'])
    for j in range(len(chunkfam.famID)):
        dosage[chunkfam.famID[j]] = pd.Series("NA", index=dosage.index)
        dosage[chunkfam.indID[j]] = pd.Series("NA", index=dosage.index)

    h.write_fam('temp/' + mf + '/'+ chunk + '_' + str(i), chunkfam)
    h.write_map('temp/' + mf + '/' + chunk + '_' + str(i), chunkmap)
    h.write_ngt('temp/' + mf + '/' + chunk + '_' + str(i), chunkngt)
    h.write_info('temp/' + mf + '/' + chunk + '_' + str(i), chunkinfo)
    dosage.to_csv('temp/' + mf + '/' + chunk + '_' + str(i) + '.gz', sep=' ', index=False, compression='gzip')



for chunk in chunks:
    wavespaths = h.read_wavespaths(s.wavesfile, chunk)

    a = []
    nwaves = len(wavespaths.filepaths)

    t = all(path.exists(wavespaths.filepaths[i] + '.gz') for i in range(nwaves))
    if not t:
        for i in range(nwaves):
            if not path.exists(wavespaths.filepaths[i] + '.gz'):
                make_nan(wavespaths, chunk, i)
                a.append('temp/' + mf + '/' + chunk + '_' + str(i) + ' ' + str(i))

    if a:
        with open('temp/' + mf + '/' + chunk + '.txt', 'w') as out:
            out.write('chunk wave\n')
            for l in a:
                out.write(l + '\n')

