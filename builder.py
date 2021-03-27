#!/usr/bin/env python

import argparse
import os
import sys
import json

sys.path.insert(0, 'modules')

parser = argparse.ArgumentParser(description="Generate DXF designs from datafiles.")

parser.add_argument('File', metavar="file-path", type=str, help="Path to the datafile file")
parser.add_argument('Output', metavar="output-path", type=str, help="Path to the output file")

args = parser.parse_args()

file_path = args.File
out_path = args.Output

try:
    dat = open(file_path, "r")
    d = dat.readlines()
except FileNotFoundError:
    print('Err: File not found! Please check your path!')
    exit()

dat = open(file_path, "r")
d2 = dat.read()

found = False

data = {}

for line in d:
    m_id = line.split()

    if m_id[0] == 'MODULE':
        module = m_id[1]
        found = True

        break

if not found:
    print('Err: Improper data file!')
    exit()

f = open("modules/tfile.dat", "a")
f.write(d2)
f.close()

os.system("python modules/nozzle.py ")
