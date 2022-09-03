#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

data = pd.read_csv(sys.argv[1], sep="\t")

data['Z'] = np.log(data.iloc[:,1])
data['Z'] = (data['Z']-data['Z'].mean())/data['Z'].std()
print(len(data))
print(np.sum([data['Z'] >= 2]))
# print(np.sum([data['Z'] >= 3]))

data[data['Z'] >= 2].to_csv(sys.argv[1][:-4] + "_Z2.tsv", sep="\t", index=False)
# data[data['Z'] >= 3].to_csv(sys.argv[1][:-4] + "_Z3.tsv", sep="\t", index=False)
