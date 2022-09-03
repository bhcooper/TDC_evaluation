#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

from pandas import read_csv

Nlookup = {"A":0, "C":1, "G":2, "T":3}

comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seqs):
    return np.array(["".join([comp[x] for x in seq][::-1]) for seq in seqs])

table = read_csv(sys.argv[1], sep="\t", index_col=0)

table_rc = table.copy()
table_rc.index = rc(table.index)
table = pd.concat((table, table_rc))
table = table.groupby(table.index).mean().sort_values(table.columns[0], ascending=False)
table.iloc[:,0] = table.iloc[:,0]/table.iloc[:,0].max()

unique = set()
for seq in table.index:
    if(not seq in unique and not rc([seq])[0] in unique):
        unique.add(seq)

table = table.loc[list(unique)]

seqs = np.array([seq.replace("N", "") for seq in table.index])
seqlen = len(seqs[0])

matrix = pd.read_csv(sys.argv[2], delim_whitespace=True, comment="#").transpose()
print(f"PWM Shape: {matrix.shape}")
matrix = np.concatenate((np.full((4, seqlen-1), 0.25), matrix, np.full((4, seqlen-1), 0.25)), axis=1)
matrix -= np.mean(matrix, axis=0).reshape(1,-1)

nshifts = matrix.shape[1] - seqlen + 1

background = 0.0**seqlen
BEEML = np.zeros((len(seqs), 2, nshifts))

for j in range(nshifts):
    BEEML[:,0,j] = [np.sum([matrix[Nlookup[c], i+j] for i,c in enumerate(seq)]) for seq in seqs]
    BEEML[:,1,j] = [np.sum([matrix[3-Nlookup[c], i+j] for i,c in enumerate(seq[::-1])]) for seq in seqs]

score = np.max(np.concatenate((np.max(BEEML, axis=2), np.full((len(seqs), 1), background)), axis=1), axis=1)
strand = np.argmax(np.concatenate((np.max(BEEML, axis=2), np.full((len(seqs), 1), background)), axis=1), axis=1)
argF = strand == 0
argR = strand == 1
argN = strand == 2
shift = np.zeros_like(strand)
shift[argF] = np.argmax(BEEML[argF,0,:], axis=1) - seqlen + 1
shift[argR] = np.argmax(BEEML[argR,1,:], axis=1) - seqlen + 1

seqs[argR] = rc(seqs[argR])
table.index = seqs
table = table.reset_index()
table.rename(columns={table.columns[0]:"Seq"}, inplace=True)

table["Shift"] = shift
table["PWM"] = score
table['background?'] = argN

minShift = -np.min(table["Shift"])
maxShift = np.max(table["Shift"])
for i, item in table.iterrows():
    shift = item["Shift"]
    table.at[i, table.columns[0]] = "N" * (minShift+shift) + item[table.columns[0]] + "N" * (maxShift-shift)

# out = table[[table.columns[0], 'PFM', table.columns[1], "Shift"]].drop_duplicates().sort_values(table.columns[1], ascending=False).reset_index(drop=True)
out = table[[table.columns[0], table.columns[1], 'PWM', "Shift"]].drop_duplicates().sort_values(table.columns[1], ascending=False).reset_index(drop=True)
out.to_csv(sys.argv[1][:-4] + "_PWM.tsv", sep="\t", index=False)
