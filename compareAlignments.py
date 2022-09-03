#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seq):
    return "".join([comp[x] for x in seq][::-1])

data1 = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
data2 = pd.read_csv(sys.argv[2], sep='\t', index_col=0)[['Shift']]

data1.index = [s.replace('N', '') for s in data1.index]
data2.index = [s.replace('N', '') for s in data2.index]

Z = data1.iloc[:,0]
Z = np.log(Z)
Z = (Z-np.mean(Z))/np.std(Z)

# print("Z = 2")
# filt = Z > 2
# data1 = data1[filt][['Shift']]

# data1[np.logical_and(filt, data1['Shift'] == 0)].to_csv(sys.argv[1][:-4] + "_Z2_shift0.tsv", sep="\t")

# data1 = data1[data1['Shift'] == 0]
# data2 = data2[data2['Shift'] == 0]

tot = len(data1)

# data1['rc'] = False
# data1_rc = data1.copy()
# data1_rc['rc'] = True
# data1_rc.index = [rc(s) for s in data1_rc.index]
# data1_rc["Shift"] = -data1_rc["Shift"]
# data1_rc = data1_rc.drop(index = data1.index.intersection(data1_rc.index))
# data1 = data1.append(data1_rc)

idxMax = data1.iloc[:,0].idxmax()
if(rc(idxMax) in data2.index):
    data2.index = [rc(s) for s in data2.index]
    data2["Shift"] = -data2["Shift"]    

data2['Shift'] += data1.loc[idxMax]['Shift'] - data2.loc[idxMax]['Shift']
data2['Shift'] = data2['Shift'].astype(int)

data = data1.merge(data2, how='inner', left_index=True, right_index=True)

print(str(np.sum(data['Shift_x'] == data['Shift_y'])) + '/' + str(tot))
exit()

if(np.sum(data['Shift_x'] == data['Shift_y']) < np.sum(data['Shift_x'] == -data['Shift_y'])):
    data['Shift_y'] = -data['Shift_y']

# output = []
# for n in range(10, 1001, 10):
# for n in range(100, len(data)+1, 100):
    # print(n)
    # output += [[n, np.sum(data['Shift_x'][:n] == data['Shift_y'][:n])/len(data[:n])]]

print(f"Matched: {np.sum(data['Shift_x'] == data['Shift_y'])}/{tot}")
# print(np.sum(np.abs(data['Shift_x']) == np.abs(data['Shift_y'])))

# output = np.array(output)
# st.saveTable("compareAlignments.tsv", output, header="# seqs\t" + sys.argv[1][:-4] + "\t" + sys.argv[2][:-4])



# matched = data['Shift_x'] == data['Shift_y']
# X = data[np.logical_or(np.logical_or(data['Shift_x'] == 0, data['Shift_y'] == 0), matched)]
# # data.to_csv(sys.argv[1].split("_")[0] + "_cores.tsv", sep="\t")

# dis = X[np.logical_not(X['Shift_x'] == X['Shift_y'])].copy()

# minShift = -np.min(dis["Shift_x"])
# maxShift = np.max(dis["Shift_x"])
# dis["Seq_x"] = ["N" * (minShift+shift) + x + "N" * (maxShift-shift) for x,shift in zip(dis.index, dis["Shift_x"])]
# minShift = -np.min(dis["Shift_y"])
# maxShift = np.max(dis["Shift_y"])
# dis["Seq_y"] = ["N" * (minShift+shift) + x + "N" * (maxShift-shift) for x,shift in zip(dis.index, dis["Shift_y"])]

# dis.to_csv("disagreements.tsv", sep="\t")