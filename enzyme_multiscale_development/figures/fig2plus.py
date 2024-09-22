#!/usr/bin/env python3
# plot histograms for a bunch of files together

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
mpl.use('svg')
import seaborn as sns
import re
import os
from cycler import cycler
#from matplotlib.ticker import FormatStrFormatter
mpl.rcParams.update({'font.size': 16})

BINCT = 50
#FILE_LIST = sys.argv[1].split()
ROOTDIR = os.getcwd()

RANGES = {
    'LIPAA': {
        'Rg':(1.9,2.2),
        'k2':(.03,.18),
        'b':(.60, 1.89)
        },
    'LIPCG': {
        'Rg':(1.8,1.95),
        'k2':(.01,.07),
        'b':(.38, 0.86)
        },
    'DHAAA': {
        'Rg':(1.7,1.9),
        'k2':(.01,.05),
        'b':(.29,.68)
        },
    'DHACG': {
        'Rg':(1.7,1.83),
        'k2':(.01,.05),
        'b':(.23, .6)
        }
}

PARAMPRINT = {
    'Rg':'$R_{g}$',
    'k2':'$\kappa^2$',
    'b':'$b$'
}

PROTPRINT = {
    'LIPAA':'Lip AA',
    'LIPCG':'Lip CG',
    'DHAAA':'Dha AA',
    'DHACG':'Dha CG',    
}

PARAMLIMITS = {
    'Rg':{'x':[1.5,2.4]},
    'k2':{'x':[0.0,0.19]},
    'b':{'x':[None,1.8]},
}

def generate_filepath(system='lipase',resolution='aa',param='Rg',seed=1):
    MYLOCATION = f'/home/mason/dmref/{system}/{resolution}/production/analysis_data/gyr/'
    MYFILE = f'{param}-seed{seed}.txt'
    return MYLOCATION+MYFILE

def make_ensemble_mean_histogram(ensemble_fileset,bin_count,hist_range,myenz,myparam,myres):
    ensemble=[]
    for myfile in ensemble_fileset:
        myfiledata = np.loadtxt(myfile, comments=('#'), usecols=(3), skiprows=1) 
        ensemble.append(myfiledata)
    np.savetxt(f'data/ensemble-{myenz}-{myparam}-{myres}.csv',ensemble,fmt='%.3f',delimiter=',')
    ensemble_flat = [ i for j in ensemble for i in j ]
    np.savetxt(f'data/ensembleflat-{myenz}-{myparam}-{myres}.csv',ensemble_flat,fmt='%.3f',delimiter=',')
    y,binEdges = np.histogram(ensemble_flat,bins=bin_count,density=True)
    if myres=='aa': ls='--'
    if myres=='cg': ls='-'
    plot_hist(y,binEdges,0,label=PROTPRINT[my_sys],linestyle=ls)

def plot_hist(counts,bins,error,label,linestyle):
    ax.stairs(counts,bins,label=label,linestyle=linestyle, lw=3)

if __name__ == '__main__':
    for param in ['Rg','k2','b']:
        for proteins in [ ['LIPAA','LIPCG', 'DHAAA','DHACG'] ]:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
            rePaired = sns.color_palette('Paired')[2:4]+sns.color_palette('Paired')[6:8]
            cc = (cycler(color=rePaired) +
                 cycler(ls=['--', '-']*2))
            ax.set_prop_cycle(cc)
            for my_sys in proteins:
                my_enz = my_sys[:3].lower()
                my_res = my_sys[3:].lower()
                if my_enz[:3] == 'lip': my_enz = 'lipase'
                if my_enz[:3] == 'dha': my_enz = 'dehalogenase'
                ensemble_fileset = [
                    generate_filepath(system=my_enz,resolution=my_res,param=param,seed=i) \
                        for i in range(1,11)
                ]
                make_ensemble_mean_histogram(ensemble_fileset,BINCT,RANGES[my_sys][param],my_enz,param,my_res)
            plt.ylim(0,)
            plt.xlim(*PARAMLIMITS[param]['x'])
            if param in ['k2']: plt.legend(frameon=False, borderpad=0)
            ax.set_xlabel(f'{PARAMPRINT[param]}')
            
#            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            plt.tight_layout()
            plt.savefig(f'images/fig2-{param}.png', dpi=300)


