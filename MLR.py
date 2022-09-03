#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd
import multiprocessing as mp

from collections import defaultdict
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle

N = ["A", "C", "G", "T"]

offset = {"A":0, "C":1, "G":2, "T":3}
def encode1mer(seq):
    X = np.zeros(len(seq)*4, dtype=bool)
    for i,c in enumerate(seq):
        if(not c == "N"): 
            X[i * 4 + offset[c]] = True
    return X

def encode1mers(seqs):
    pool = mp.Pool(processes = mp.cpu_count())
    X = pool.map(encode1mer, seqs)
    pool.terminate()
    return np.array(X)

def encodeShape(seqs, shape):
    data = pd.read_csv(os.path.dirname(sys.argv[0]) + "/featureTables/" + shape + ".tsv", sep='\t', usecols=(0,1))
    kmers = data.iloc[:,0].values
    y = data.iloc[:,1].values
    featLen = len(kmers[0])
    lookup = defaultdict(lambda:np.mean(y), zip(kmers, y))
    X = [[lookup[seq[i:i+featLen]] for i in range(len(seq)-featLen+1)] for seq in seqs]
    X = np.array(X)
    return X

data = pd.read_csv(sys.argv[1], sep="\t", usecols=(0,1))
seqs = data.iloc[:,0].values
y = data.iloc[:,1].values 

argsort = np.argsort(-y)
seqs = seqs[argsort]
y = y[argsort]
seqs, y = shuffle(seqs, y, random_state=0)

y = np.log(y)

X_1mer = encode1mers(seqs)
X_mgw = encodeShape(["NN" + s + "NN" for s in seqs], "MGW")
X_EP = encodeShape(["NN" + s + "NN" for s in seqs], "EP")

X = np.concatenate((X_1mer, X_mgw, X_EP), axis=1)
X = StandardScaler().fit_transform(X)

model = ElasticNetCV(max_iter = 1e6, cv=5, n_jobs=-1)
model = cross_validate(model, X, y, cv=5, scoring="r2", return_train_score=True)
print('Train R²: %.4f' % np.median(model["train_score"]))
print('Test  R²: %.4f' % np.median(model["test_score"]))
