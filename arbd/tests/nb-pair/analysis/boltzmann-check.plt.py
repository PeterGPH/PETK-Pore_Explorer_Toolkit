#! /usr/bin/env python3
# /home/cmaffeo2/bin/python3.sh
# /software/python3-3.4.1/bin/python3

import numpy as np
from scipy import exp
from scipy import sqrt
import matplotlib.pyplot as plt
# import scipy.optimize as opt
from glob import glob
# from natsort import natsorted
import re

from plotTools import *

import sys
def info(*obj):
    print('INFO: ',*obj , file=sys.stderr)

def makeTitle(ax,text):
    # ax.set_title(text)
    ax.annotate(text, xy=(0.05,0.92), 
                fontsize=9, xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top')
    
fig,axes = plt.subplots(1)
ax1 = axes

makeTitle(ax1,'Pair interaction (rev. a158e71)')
ax1.set_ylabel('likelihood')
ax1.get_yaxis().set_label_coords(-0.085,0.5)
ax1.set_xlabel("distance between isolated pair of particles")

xmin, xmax = [2,8]
setRange(ax1,[xmin,xmax, 0,0.2])

def loadData(fname):
    d0 = np.loadtxt(fname)
    x0 = d0[:,0]
    y0 = d0[:,1]
    
    ids = np.where( (x0 > xmin) & (x0 < xmax) )
    x0 = x0[ids]
    y0 = y0[ids]

    y0 = y0 / np.sum(y0)

    return x0,y0

kT=0.58622592
def estimate():
    eps = 2*10
    sig = 3

    x = np.linspace(xmin,xmax,100)
    u = 4*eps*( (sig/x)**12 - (sig/x)**6 )
    
    y = np.exp( - u / kT ) * x**2
    y = y / np.sum(y)
    return x,y

def estimate():
    d0 = np.loadtxt('../dummy-pot.dat')
    x = d0[10:,0]
    u = d0[10:,1]

    ids = np.where( (x > xmin) & (x < xmax) )
    x = x[ids]
    u = u[ids]

    y = np.exp( - u/kT ) * x**2

    y = y / np.sum(y)
    return x,y


def plot(x,y,l,c):
    # return ax1.plot(x,y, ls='-', marker='o', color=c, mec=c)
    return ax1.bar(x,y, label=l, lw=0, width=np.mean(x[1:]-x[:-1]), alpha=0.5, color=c)


x,y = loadData('rho.dat')
ex,ey = estimate()

plot(x,y,'simulated', c0)
plot(ex,ey,'expected',c1)

plt.legend()

# plt.show()
plt.tight_layout()
# plt.savefig(sys.stdout, format='svg')
plt.savefig(sys.argv[1])
