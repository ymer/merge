import pandas as pd
from glob import glob
import helper as h
from os import path


s = h.parse()


def collect_tests():
    with open(s.checkfolder + "check.txt", 'w') as out:
        paths = glob(s.checkfolder + "chunks/*_check.txt")
        for pth in paths:
            df = pd.read_csv(pth, sep="\s+", header=None)
            if any(df[0] != True):
                out.write(path.basename(pth) + "\n")


def collect_missing_stats():
    paths = glob(s.checkfolder + "chunks/*_missing_stats.txt")

    df = pd.read_csv(paths[0], sep="\t")
    for i in range(1, len(paths)):
        ndf = pd.read_csv(paths[i], sep="\t")
        df = df.add(ndf)

        df["present_in_waves"] = ndf["present_in_waves"]

    df.to_csv(s.checkfolder + "missing_stats.txt", sep="\t", index=False)


def collect_compare():
    with open(s.checkfolder + "compare.txt", 'w') as out:
        paths = glob(s.checkfolder + "compare/*.txt")
        for path in paths:
            df = pd.read_csv(path, sep="\t", header=None)
            if any(df[0] != True):
                out.write(path + "\n")


def collect_missing_with_indco():

    files = glob(s.checkfolder + "chunks/*_missing_with_indco.txt")

    with open(s.checkfolder + "missing_with_indco.txt","w") as outfile:
        with open(files[0]) as f1:
            for line in f1:
                outfile.write(line)

        for x in files[1:]:
            with open(x) as f1:
                for line in f1:
                    if not line.startswith("SNP"):
                        outfile.write(line)


def collect_missing():

    files = glob(s.checkfolder + "chunks/*_missing.txt")

    with open(s.checkfolder + "missing.txt","w") as outfile:
        with open(files[0]) as f1:
            for line in f1:
                outfile.write(line)

        for x in files[1:]:
            with open(x) as f1:
                for line in f1:
                    if not line.startswith("SNP"):
                        outfile.write(line)


def collect_numbers():

    files = glob(s.checkfolder + "chunks/*_missing.txt")

    totals = []
    for f in files:
        df = pd.read_csv(f, sep='\s+')
        total = df.sum(numeric_only=True)
        total.name = path.basename(f)[:-12]
        totals.append(total)

    snp_totals = pd.concat(totals, axis=1).transpose()

    snp_totals.to_csv(s.checkfolder + 'snp_totals.txt', sep=' ')


collect_tests()

collect_missing_stats()

collect_missing()

collect_numbers()

collect_missing_with_indco()

