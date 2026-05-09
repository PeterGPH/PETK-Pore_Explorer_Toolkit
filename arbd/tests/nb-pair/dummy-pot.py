#! /usr/bin/env python

import numpy as np

eps=3
sig=3

r = np.linspace(0,50,501)
r[0] = 0.1
u = 4*eps*( (sig/r)**12 - (sig/r)**6 )
r[0] = 0

ids = np.where( u > 1000 )
u[ids] = u[ids[0][-1]]

u = u + r;
u = u-u[-1]

ch = open('dummy-pot.dat','w')

for data in zip(r,u):
    ch.write("%f %f\n" % data)
