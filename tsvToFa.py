#!/usr/bin/env python

import sys

outfile = open(sys.argv[1][:-4] + ".fa", 'w')
pad = "N" * 9

for i,line in enumerate(open(sys.argv[1])):
    if i > 0:
        outfile.write(">"+str(i)+"\n" + pad + line.split()[0].replace("_", "") + pad + "\n")

outfile.close()