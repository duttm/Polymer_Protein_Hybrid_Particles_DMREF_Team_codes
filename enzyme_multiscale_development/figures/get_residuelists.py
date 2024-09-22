#!/usr/bin/env python3

# get a list of surface residues as identified in surfmap output
# to be used in preparation of molecular visualizations

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
mpl.use('svg')
import seaborn as sns
import re
import os
from cycler import cycler

ROOTDIR = os.getcwd()

rescounts = {}
p=re.compile(r'[A-Z]+_.*_A')
myfile = sys.argv[1]
with open(myfile) as f:
    for line in f:
        matches = re.findall(p,line)
        if matches: 
            splitmatches = matches[0].split(',')
#            print(splitmatches)
            for i in splitmatches:
                matchdata = i.strip().split('_')
#                print(matchdata)
                resname = matchdata[0]
                resnum = int(matchdata[1])
#                print(resnum)
                if resname not in rescounts.keys():
                    rescounts[resname] = set()
#                    print(rescounts)
                rescounts[resname].update([resnum])
    #        print(rescounts)

for i in sorted(rescounts.keys()):
    print(i, ' '.join(str(mynum) for mynum in sorted(rescounts[i])))

#print(rescounts)
