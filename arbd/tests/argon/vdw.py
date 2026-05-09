#! /usr/bin/env python

import numpy as np

eps=0.238000
rmin = 2*1.908100
sig = rmin/(2**(1/6))

r = np.linspace(0,50,501)
r[0] = 0.1
u = 4*eps*( (sig/r)**12 - (sig/r)**6 )
r[0] = 0

ids = np.where( u > 200 )
idMax = ids[0][-1]
du = u[idMax] - u[idMax-1]
u[ids] = u[idMax] -du * (np.arange(idMax,-1,-1)  )

u = u-u[-1]

ch = open('dummy-pot.dat','w')
for data in zip(r,u):
    ch.write("%f %f\n" % data)
