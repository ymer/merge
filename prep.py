import pandas as pd

bn = '/project/DBS1to23/imatthei_imputation/phase3_woBLACKLIST_newQC/'
pn1 = 'dasuqc1_dbs_'
pn2 = '_mix_mm-qc.hg19.ch.fl/'
fn1 = 'dos_dbs_dbs'
fn2 = '_mix_mm-qc.hg19.ch.fl.'
w = 'dbs'
nwaves = 23
filename = 'files/file_locations_phase3_woBLACKLIST_newQC.txt'

d = {}
for wave in range(1, nwaves + 1):
    wave = str(wave)
    path = bn + w + wave + '/' + pn1 + w + wave + pn2 + 'qc1/' + fn1 + wave + fn2
    outname = fn1 + 'ALL' + fn2
    infopath = bn + w + wave + '/' + pn1 + w + wave + pn2 + 'info/' + fn1 + wave + fn2
    wavename = 'w' + wave
    wavefolder = pn1 + w + wave + pn2
    d[int(wave) -1] = {'infopath': infopath, 'outname': outname, 'path': path, 'wavefolder': wavefolder, 'wavename': wavename}

df = pd.DataFrame(d).transpose()

df.to_csv(filename, sep=' ')




