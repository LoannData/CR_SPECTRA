#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 08:29:50 2020

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as sc

import sys 
sys.path.append("../")
import constants as cst
import phases_collection as ism 

phase = ism.CNM
B0 = phase.get('B')
mi = phase.get("mi")
mn = phase.get("mn")
ni = phase.get("ni")
nn = phase.get("nn")
T  = phase.get("T")

chi = (mn*nn)/(mi*ni)
VAi = B0/np.sqrt(4*np.pi*mi*ni)
VA    = B0/np.sqrt(4*np.pi*(mi*ni + mn*nn))
if (T <= 50) : 
    nu_in = 2*nn*8.4e-9*(50/1e4)**0.4
if (T > 50) : 
    nu_in = 2*nn*8.4e-9*(T/1e4)**0.4
nu_ni = chi**(-1.)*nu_in


k_min = 1e-20 # [cm]
k_cm  = 1e-15 # [cm]
k_cp  = 1e-20 # [cm]
k_max = 2*nu_ni/VA  # [cm] 


def J(a, b, c, D, q, x) : 
    d = 2*a*x + b
    A1 = 1./(q*D)*2**(q)*x**(-q)
    A2 = (a*x)/(D + d)**q
    A3 = sc.hyp2f1(q, q, q+1, (b + D)/(d + D))
    A4 = (a*x)/(- D + d)**q
    A5 = sc.hyp2f1(q, q, q+1, (b - D)/(d - D))
    return A1*(A2*A3 - A4*A5)



E = 10*cst.GeV
m = cst.mp
gamma = 1 + (E /(m*cst.c**2))
v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
p = gamma*m*v 
Omega0 = cst.e*B0/(m*cst.c)
Omega  = Omega0/gamma 

mu = np.linspace(-0.99, 0.99, 100)
d_uu  = np.zeros(len(mu))
k_zz = 0. 
dmu = mu[1] - mu[0]

Itot = 1e-1#(1e-2)**2*(k_max - k_cp)

for ii in range(len(mu)) : 
    q = 5./3

    a = (v*mu[ii] - VA)**2
    b = 2*Omega*(v*mu[ii] - VA)
    c = Omega**2 + (- nu_in/2)**2
    #d = 2*a*mu[ii] + b 
    D = b**2 - 4*a*c
    
    eps = VA/v
    gs0 = 2*(q-1)*(B0**2/(8*np.pi))*Itot*k_cp**(q - 1)
    
    
    
    
    fj_p = J(a, b, c,  D, q, k_max) - J(a, b, c,  D, q, k_cp)
    fj_m = J(a, -b, c,  D, q, k_max) - J(a, -b, c,  D, q, k_cp)
    
    I = -2*gs0*(1 - mu[ii]*eps)**2*(- nu_in/2.)*fj_p
    
    d_uu[ii] = np.pi*Omega**2*(1 - mu[ii]**2)/B0**2*I
    
    if (d_uu[ii] > 0.) : 
        k_zz += v**2/8.*(1 - mu[ii]**2)**2/d_uu[ii]*dmu




plt.figure()
plt.semilogy(mu, d_uu)



print ("k_zz = ",k_zz," cm^2/s")
























