#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 15:31:03 2019

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 


GeV     = 0.00160218
kyr     = 1e3*24*60*60*365.25

rho_0   = 0.35#1.24e-24   # [1/cm^3]
e       = 4.8032e-10 # [statcoul] 
c       = 2.998e10   # [cm/s]
xhi_cr  = 0.1        # CR acceleration efficiency 
xhi_0   = 2.026
beta    = 0.2 
Esn     = 1e51       # Energy liberated by the SN
Emin    = 0.1*GeV   # Minimum energy produced during sedov phase 

def tesc(E) : 
    A = 4*np.sqrt(np.pi*rho_0)*e/(125*c)*xhi_cr
    B = (xhi_0*Esn/rho_0)**(3./5)
    C = beta/(1 + beta)/E
    D = Emin**beta/(E**beta - Emin**beta)
    return (A*B*C*D)**(5./4)

def gauss(t, sig, mu) : 
    return np.exp(-0.5*((t - mu)/sig)**2)

def sig(t) : 
    return t/10.


E = np.logspace(np.log10(10.*GeV), np.log10(1e3*GeV), 100)
#t = np.logspace(np.log10(0.01*kyr), np.log10(1e3*kyr), 1000)
t = np.linspace(0.01*kyr, 1e3*kyr, 10000)
#sig = 0.1*kyr


Qcr = np.empty((len(E), len(t)))
for ei in range(len(E)) : 
    for ti in range(len(t)) : 
        Qcr[ei][ti] = gauss(t[ti], sig(t[ti]), tesc(E[ei]))
    Qcr[ei] = Qcr[ei]#/sum(Qcr[ei])

plt.figure(figsize=(10,6))
plt.semilogx(t/kyr, Qcr[0], label="E = "+str(round(E[0]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[10], label="E = "+str(round(E[10]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[20], label="E = "+str(round(E[20]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[30], label="E = "+str(round(E[30]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[40], label="E = "+str(round(E[40]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[50], label="E = "+str(round(E[50]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[60], label="E = "+str(round(E[60]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[70], label="E = "+str(round(E[70]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[80], label="E = "+str(round(E[80]/GeV,1))+" GeV")
plt.semilogx(t/kyr, Qcr[90], label="E = "+str(round(E[90]/GeV,1))+" GeV")

plt.legend()
plt.xlabel("Time [kyr]")
plt.ylabel("Injection rate")

plt.savefig("Injection.pdf")


#t = np.empty(len(E))
#
#for ii in range(len(E)) :
#    t[ii] = tesc(E[ii])
#
#
#plt.figure(figsize=(8,6))
#plt.loglog(E/GeV, t/kyr)
