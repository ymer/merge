import pandas as pd
import helper as h

s = h.parse()
s.wavepaths = h.replace(s)

wave_fams = h.read_wave_fams(s.wavepaths, usecols=["famID", "indID"])

for i, fam in enumerate(wave_fams):
    fam["wave"] = i + 1

merged_fam = h.read_fam(s.mergepath, usecols=["famID", "indID"])

waves = pd.concat(wave_fams)

merged_fam = pd.merge(merged_fam, waves, how='left', on=['famID', 'indID'])

merged_fam.to_csv(s.checkfolder + "covar.txt", sep=" ", header=False, index=False)


