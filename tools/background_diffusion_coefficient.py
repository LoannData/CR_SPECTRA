#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 09:17:25 2020

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as sc
import math 

# import sys 
# sys.path.append("../")
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


def Duu_Alfven_Slab_Linear_Undamped(mu, E, medium_props, mass = cst.mp, kmin = 1e-20, q = 1.5, I = 1e-4) : 
    """
    See. Schlickeiser (2002, Chap 13.1.3.1, p.318)

    Parameters
    ----------
    mu : TYPE
        DESCRIPTION.
    medium_props : TYPE
        DESCRIPTION.
    particles_props : TYPE
        DESCRIPTION.
    kmin : TYPE, optional
        DESCRIPTION. The default is 1e-20.
    q : TYPE, optional
        DESCRIPTION. The default is 1.5.

    Returns
    -------
    D : TYPE
        DESCRIPTION.

    """
    
    def g_s(k, kmin, q, medium_props, I = I) : 
        B0 = medium_props.get("B")
        W0 = B0**2/(8*np.pi)
        I  = I
        g_s0 = 2*(q - 1)*W0*I*kmin
        if (k >= 0) : 
            if (k >= kmin) : 
                return g_s0*k**(-q)
            else : 
                return 0. 
        if (k < 0) : 
            if (abs(k) >= kmin) : 
                return g_s0*abs(k)**(-q)
            else : 
                return 0.
    
    
    m = mass 
    gamma = 1 + (E /(m*cst.c**2))
    v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
    p = gamma*m*v 
    Omega0 = cst.e*medium_props.get("B")/(m*cst.c)
    Omega  = Omega0/gamma 
    
    
    B0    = medium_props.get("B")
    mn    = medium_props.get("mn")
    nn    = medium_props.get("nn")
    mi    = medium_props.get("mi")
    ni    = medium_props.get("ni")
    
    # eps = B0/np.sqrt(4*np.pi*(mi*ni + mn*nn))/v
    eps = B0/np.sqrt(4*np.pi*(mi*ni))/v
    
    k_rp = (Omega/v)/(mu - eps)
    k_rm = (Omega/v)/(mu + eps)
    
    Cte1    = np.pi*Omega**2*(1 - mu**2)/(v*B0**2)
    Cte2_jp = (1 - eps*mu)**2/abs(mu - eps)
    Cte2_jm = (1 + eps*mu)**2/abs(mu + eps)
    Gp_jp      = g_s(k_rp, kmin, q, medium_props)  
    Gm_jp      = g_s(-k_rp, kmin, q, medium_props) 
    Gp_jm      = g_s(k_rm, kmin, q, medium_props)  
    Gm_jm      = g_s(-k_rm, kmin, q, medium_props) 
    
    Djp = Cte1*Cte2_jp*(Gp_jp + Gm_jp)
    Djm = Cte1*Cte2_jm*(Gp_jm + Gm_jm)
    
    D = Djp + Djm
    return D




def Kappa_zz(E, medium_props, mass = cst.mp, kmin = 1e-20, q = 5./3, I = 1e28) : 
    """
    

    Parameters
    ----------
    E : TYPE
        DESCRIPTION.
    medium_props : TYPE
        DESCRIPTION.
    mass : TYPE float, optional
        Mass of the diffusin particle. The defaut is m_proton
    kmin : TYPE float, optional
        Minimun length in cm^-1 for the turbulence spectra. The default is 1e-20.
    q : TYPE float, optional
        Spectral index of the Kolmogorov-like turbulence 
        spectrum. The default is 5./3.
    I : TYPE float, optional
        Diffusion coefficient normalization value for 1GeV particle and 5./3 spectrum. The default is 1e28.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    I = (2.2345649450772952e28*1e6/I)
    
    m = mass
    gamma = 1 + (E /(m*cst.c**2))
    v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
    p = gamma*m*v 
    Omega0 = cst.e*medium_props.get("B")/(m*cst.c)
    Omega  = Omega0/gamma 
    
    k_zz = 0.
    mu = np.linspace(-1, 1, 100)
    for jj in range(1, len(mu)-1) : 
        dmumu = 0.5*(mu[jj+1] - mu[jj-1])
        k_zz += dmumu*(1 - mu[jj]**2)**2/Duu_Alfven_Slab_Linear_Undamped(mu[jj], E, medium_props, mass = m, kmin = kmin, q = q, I = I)
    
    return v**2/8*k_zz




# Experiment 
# Emin = 0.1*cst.GeV
# Emax = 100*cst.TeV
# E = np.logspace(np.log10(Emin), np.log10(Emax), 100)
# K_zz_1 = np.zeros(len(E))
# K_zz_2 = np.zeros(len(E))

# for ii in range(len(E)) : 
#     K_zz_1[ii] = Kappa_zz(E[ii], ism.WNM, mass = cst.mp, kmin = 1e-20, q = 5./3, I = 1e28)
#     K_zz_2[ii] = Kappa_zz(E[ii], ism.WNM, mass = cst.mp, kmin = 1e-20, q = 5./3, I = 1e28)

# plt.loglog(E/cst.GeV, K_zz_1, c="green")
# plt.loglog(E/cst.GeV, K_zz_2, c="blue")

