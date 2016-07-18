from subprocess import call
import argparse
from sys import argv
import pickle


parser = argparse.ArgumentParser()

parser.add_argument('--gwf', action='store', dest='gwf', help='path to the gwf program')
parser.add_argument('--date', action='store', dest='date', help='data (or other identifier)')
parser.add_argument('--vers', action='store', dest='vers', help='version of the data')
parser.add_argument('--chunks', action='store', dest='chunks', default='one', help='The number of chunks to merge (all / sample/ one) (default: one)')
parser.add_argument('--genoprob', action='store', dest='genoprob', default=0.0, type=float, help='genoprob cutoff')
parser.add_argument('--run_check', action='store_true', dest='run_check', default=True, help='If set, perform a check of the merging')
parser.add_argument('--indlist', action='store', dest='indlist', default=False, help='File with list of individuals to include (default: include all)')

s = parser.parse_args(argv[1:])

pickle.dump(s, open('files/config.pl', 'wb'))

call([s.gwf])


