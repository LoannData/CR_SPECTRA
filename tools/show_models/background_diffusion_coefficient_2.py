#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 12:50:39 2020

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as sc
import math 

import sys 
sys.path.append("../")
import constants as cst
import phases_collection as ism 
import mathmethods as mt 









def IonNeutral_Damping(k, medium_props, nu_n = 0, theta = 0) : 
    T  = medium_props.get('T')
    B  = medium_props.get('B')
    mn = medium_props.get('mn')
    mi = medium_props.get('mi') 
    X  = medium_props.get('X')
    ni = medium_props.get('ni')
    nn = medium_props.get('nn')    
    
    
    # Some relations 
    # rl = (E)/(cst.e*B)
    # k = 1/rl
    chi = (mn/mi)*(X**(-1.)-1.)
    
    VAi = B/np.sqrt(4*np.pi*mi*ni)
    if (T <= 50) : 
        nu_in = 2*nn*8.4e-9*(50/1e4)**0.4
    if (T > 50) : 
        nu_in = 2*nn*8.4e-9*(T/1e4)**0.4
    nu_ni = chi**(-1.)*nu_in
    
    kz = k*np.cos(theta)
    
    a = 1 
    b = (1 + chi)*nu_ni
    c = kz**2*VAi**2
    d = nu_ni*kz**2*VAi**2
    
    roots = mt.Cubic3(a, b, c, d)
    
    wR = abs(roots[2].imag)
    wI = abs(roots[2].real)
    cA = wR/k 
    
    return dict(wr=wR, wi=-wI, VA=cA)

medium_props = ism.WNM


E = 10*cst.GeV
m = cst.mp
gamma = 1 + (E /(m*cst.c**2))
v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
p = gamma*m*v 
Omega0 = cst.e*medium_props.get("B")/(m*cst.c)
Omega  = Omega0/gamma 

particles_props = {"v":v, "Omega":Omega}

def R(k, medium_props, particles_props, mu, kind = "+") : 
    w = IonNeutral_Damping(k, medium_props)
    v     = particles_props.get("v")
    Omega = particles_props.get("Omega")
    
    if (kind == "+") : 
        return w.get("wi")/(w.get("wi")**2 + (v*mu*k - w.get("wr") + Omega)**2)
    if (kind == "-") : 
        return w.get("wi")/(w.get("wi")**2 + (v*mu*k - w.get("wr") - Omega)**2)



def g_s(k, kmin, q, medium_props) : 
    B0 = medium_props.get("B")
    W0 = B0**2/(8*np.pi)
    I  = 1e-4
    g_s0 = 2*(q - 1)*W0*I*kmin
    if (k >= kmin) : 
        return g_s0*k**(-q)
    else : 
        return 0. 


mu = np.linspace(-0.99, 0.99, 100)
k = np.logspace(-20, -10, 100)
kmin = 1e-17
q = 1.5


k_zz = 0.

for jj in range(len(mu)) : 
    D_uu = 0. 
    for ii in range(1, len(k)-1) : 
        dk = 0.5*(k[ii+1] - k[ii-1])
        
        Omega = particles_props.get("Omega")
        B0    =  medium_props.get("B")
        
        w = IonNeutral_Damping(k[ii], medium_props)
        wtot = w.get("wr") + 1j*w.get("wi")
        
        I = dk*g_s(k[ii], kmin, q, medium_props)*(1 - (mu[jj]*wtot)/(k*v))**2
        I = I*(R(k[ii], medium_props, particles_props, mu[jj], kind = "+") + R(k[ii], medium_props, particles_props, mu[jj], kind = "-"))
        
        D_uu += np.pi*Omega**2*(1 - mu[jj]**2)/B0**2
    k_zz += v**2/8*(1 - mu[jj]**2)**2/D_uu
    
    
print ("k_zz = ",np.log10(k_zz)," cm^2/s")

# k = np.logspace(-20, -10, 1000)
# roots = []
# wR = np.zeros(len(k))
# wI = np.zeros(len(k))
# cA = np.zeros(len(k))


# for ii in range(len(k)) : 
#     w =  IonNeutral_Damping(k[ii] ,ism.WNM, nu_n = 0, theta = 0)
#     wR[ii] = w.get("wr")
#     wI[ii] = w.get("wi")
#     cA[ii] = w.get("VA")


# plt.figure()
# plt.loglog(k, wR)
# plt.loglog(k, -wI)
# # plt.loglog(k, cA)


