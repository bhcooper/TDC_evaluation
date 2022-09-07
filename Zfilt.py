#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from collections import defaultdict

comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seqs):
    return np.array(["".join([comp[x] for x in seq][::-1]) for seq in seqs])

def consolidate(lookup):
    for key in list(lookup.keys()):
        if(lookup[rc([key])[0]] == 0.0):
            lookup[rc([key])[0]] = lookup[key]
        else:
            ave = (lookup[key] + lookup[rc([key])[0]]) / 2
            lookup[key] = ave
            lookup[rc([key])[0]] = ave
    return lookup

data = pd.read_csv(sys.argv[1], sep="\t")
seqs = data.iloc[:,0]
scores = data.iloc[:,1]
lookup = defaultdict(float, zip(seqs, scores))
lookup = consolidate(lookup)
data.iloc[:,1] = [lookup[seq] for seq in seqs]


unique = set()
for seq in data.iloc[:,0]:
    if(not rc([seq])[0] in unique):
        unique.add(seq)

data = data.set_index(data.columns[0], drop = True)
data = data.loc[list(unique)].sort_values(data.columns[0], ascending=False).copy()

data['Z'] = np.log(data.iloc[:,0])
data['Z'] = (data['Z']-data['Z'].mean())/data['Z'].std()

data[data['Z'] >= 2].to_csv(sys.argv[1][:-4] + "_Z2.tsv", sep="\t", index=True)
