#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
import os
#mpl.use('svg') # svg backend for output to svg
mpl.use('GTK4Agg') # gui backend for debugging
import seaborn as sns
np.set_printoptions(threshold=sys.maxsize)
from pprint import pprint
from collections import Counter
import scipy.stats
from cycler import cycler
import matplotlib.patches as mpatches
from pprint import pprint
mpl.rcParams.update({'font.size': 15})



ROOTDIR = os.getcwd()
DATADIR='/path/to/data/'
#MYFILE=sys.argv[1]

data = [
    'files.txt',
    ]

def load_and_norm(filename):
    mydata = np.loadtxt(filename, comments=('@','#'))
    bounds = (mydata[mydata[:,1]!=0][0][0], mydata[mydata[:,1]!=0][-1][0])
    simdatanorm = mydata[:,1] / np.max(mydata[:,1])
#    print(mydata)
#    print(mydata[:,0])
#    print(mydata[:,1])
    normdata = np.column_stack((mydata[:,0],simdatanorm))
    maxloc = mydata[mydata.argmax(axis=0)[1]][0]
    return normdata, bounds, maxloc

rePaired = sns.color_palette('Paired')[0:2]+sns.color_palette('Paired')[4:6]+sns.color_palette('Paired')[8:10]
custom_cycler = (cycler(color=rePaired) +
                 cycler(ls=['dashed', 'solid']*2+['solid','solid']))

fig,ax = plt.subplots(1,1,figsize=(6,9/16*6))
ax.set_prop_cycle(custom_cycler)
xlim = [None,None]
for file in data:
    filepath = DATADIR + file
    normdata, bounds, maxloc = load_and_norm(filepath)
    print(maxloc)
    if xlim[0]:
        if bounds[0] < xlim[0]: xlim[0] = bounds[0]
    if xlim[1]: 
        if bounds[1] < xlim[1]: xlim[1] = bounds[1]
    ax.plot(normdata[:,0],normdata[:,1],lw=3)

ax.set_xlim(xlim)
#ax.set_xlim([0,5])
#ax.set_ylim([None,0.2])
plt.legend(['Lip SAXS','Lip AA','Dha SAXS','Dha AA'])
plt.xlabel('$r$ [$\AA$]')
plt.ylabel('$P(r)/P_{max}(r)$')
#plt.title(f'Pair-Distance Distribution Function P(r)\nsimulation vs experiment')
plt.tight_layout()
#plt.show()
plt.savefig('images/fig3.png')


