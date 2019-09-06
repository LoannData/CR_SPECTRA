#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 22:44:55 2019

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 

def gauss(t, sig, mu) : 
    return np.exp(-0.5*((t-mu)/sig)**2)




tmin = 0.01
tesc = 25.23
tmax = 2*tesc - tmin 

sig = 2
mu = 25

# Code resolution
Nc = 1e6
tc = np.linspace(tmin, tmax, Nc)

# Approx resolution 
Nv = 100
ta = np.linspace(tmin, tmax, Nv)

# Ratio
r = Nc/Nv 

gauss_c = np.zeros(len(tc))
gauss_a = np.zeros(len(ta))

C_c = 0. 
C_a = 0. 

for ti in range(len(tc)) : 
    C_c += gauss(tc[ti], sig, mu)

for ti in range(len(ta)) : 
    C_a += gauss(ta[ti], sig, mu)

for ti in range(len(tc)) : 
    gauss_c[ti] = gauss(tc[ti], sig, mu)/(C_a*r)

for ti in range(len(ta)) : 
    gauss_a[ti] = gauss(ta[ti], sig, mu)/C_a



plt.plot(tc, gauss_c, color = "red")
plt.plot(ta, gauss_a, color = "blue")
plt.axvline(tesc, c="black")
plt.legend()