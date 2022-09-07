#!/usr/bin/env python3

import sys
import gzip
import numpy as np

# Read sequence lines from gzip file
seqs = np.array([line.strip() for i, line in enumerate(gzip.open(sys.argv[1], "rt", encoding='utf-8')) if i % 4 == 1])

unique,counts = np.unique(seqs, return_counts=True)

# Combine seqs and counts into a 2D array
temp = np.array([unique, counts]).transpose()
# Sort by counts, descending
temp = temp[np.argsort(-counts)]

np.savetxt(sys.argv[1][:-9] + ".tsv", temp, delimiter="\t", header="", comments="", fmt="%s")