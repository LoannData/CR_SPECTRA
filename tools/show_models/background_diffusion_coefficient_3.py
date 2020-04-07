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

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys 
sys.path.append("../")
import constants as cst
import phases_collection as ism 
import mathmethods as mt 
import damping as dp 


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
        g_s0 = 2*(q - 1)*W0*I*kmin**(q - 1)
        # g_s0 = I*W0*kmin
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
    
    # I = (2.2345649450772952e28*1e6/I)
    
    
    
    m = mass
    gamma = 1 + (E /(m*cst.c**2))
    v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
    p = gamma*m*v 
    Omega0 = cst.e*medium_props.get("B")/(m*cst.c)
    Omega  = Omega0/gamma 
    
    lz_damping = dp.damping_lazarian_nopos(E, medium_props)
    # print (lz_damping)
    kmax = (lz_damping[1])**(-1)
    # kmax = 1e-21
    
    # print (kmax)
    
    I = I*kmax
    
    k_zz = 0.
    mu = np.linspace(-1, 1, 100)
    for jj in range(1, len(mu)-1) : 
        dmumu = 0.5*(mu[jj+1] - mu[jj-1])
        k_zz += dmumu*(1 - mu[jj]**2)**2/Duu_Alfven_Slab_Linear_Undamped(mu[jj], E, medium_props, mass = m, kmin = kmin, q = q, I = I)
    
    return v**2/8*k_zz



def kappa_zz_BC(E, d00, delta) : 
    return d00*(E/(10*cst.GeV))**(delta)







# Experiment 
Emin = 0.1*cst.GeV
Emax = 100*cst.TeV
E = np.logspace(np.log10(Emin), np.log10(Emax), 100)

K_zz_bc_1 = np.zeros(len(E))
K_zz_bc_2 = np.zeros(len(E))

K_zz_1    = np.zeros(len(E))
K_zz_2    = np.zeros(len(E))

K_zz_HII_1 = np.zeros(len(E))
K_zz_WIM_1 = np.zeros(len(E))
K_zz_WNM_1 = np.zeros(len(E))
K_zz_CNM_1 = np.zeros(len(E))
K_zz_DiM_1 = np.zeros(len(E))
K_zz_DeM_1 = np.zeros(len(E))
K_zz_DeC_1 = np.zeros(len(E))

K_zz_HII_2 = np.zeros(len(E))
K_zz_WIM_2 = np.zeros(len(E))
K_zz_WNM_2 = np.zeros(len(E))
K_zz_CNM_2 = np.zeros(len(E))
K_zz_DiM_2 = np.zeros(len(E))
K_zz_DeM_2 = np.zeros(len(E))
K_zz_DeC_2 = np.zeros(len(E))

for ii in range(len(E)) : 
    K_zz_bc_1[ii] =  kappa_zz_BC(E[ii], 1e28, 0.5)
    K_zz_bc_2[ii] =  kappa_zz_BC(E[ii], 1e29, 2./3)
    
    # K_zz_1[ii] = Kappa_zz(E[ii], ism.WNM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    # K_zz_2[ii] = Kappa_zz(E[ii], ism.HII, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    
    K_zz_HII_1[ii] = Kappa_zz(E[ii], ism.HII, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    K_zz_WIM_1[ii] = Kappa_zz(E[ii], ism.WIM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    K_zz_WNM_1[ii] = Kappa_zz(E[ii], ism.WNM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    K_zz_CNM_1[ii] = Kappa_zz(E[ii], ism.CNM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    K_zz_DiM_1[ii] = Kappa_zz(E[ii], ism.DiM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    K_zz_DeM_1[ii] = Kappa_zz(E[ii], ism.DeM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    K_zz_DeC_1[ii] = Kappa_zz(E[ii], ism.DeC, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 1.5, I = 1e-2)
    
    K_zz_HII_2[ii] = Kappa_zz(E[ii], ism.HII, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    K_zz_WIM_2[ii] = Kappa_zz(E[ii], ism.WIM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    K_zz_WNM_2[ii] = Kappa_zz(E[ii], ism.WNM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    K_zz_CNM_2[ii] = Kappa_zz(E[ii], ism.CNM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    K_zz_DiM_2[ii] = Kappa_zz(E[ii], ism.DiM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    K_zz_DeM_2[ii] = Kappa_zz(E[ii], ism.DeM, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    K_zz_DeC_2[ii] = Kappa_zz(E[ii], ism.DeC, mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)




size_x = 4
size_y = 3
sub_x  = 2
sub_y  = 2
fig = plt.figure(figsize=(size_x*sub_x,size_y*sub_y))

gs = gridspec.GridSpec(ncols= sub_x, nrows = sub_y, figure = fig )
gs.update(wspace=0.05, hspace=0.05) # set the spacing between axes.


ax0 = fig.add_subplot(gs[0, 0])

ax0.loglog(E/cst.GeV, K_zz_bc_1, c="black", ls = '--')
ax0.loglog(E/cst.GeV, K_zz_bc_2, c="black")
ax0.fill_between(E/cst.GeV, K_zz_bc_1, K_zz_bc_2, facecolor = "black", alpha = 0., hatch="\\")

ax0.loglog(E/cst.GeV, K_zz_HII_1, c="red", ls = "--")
ax0.loglog(E/cst.GeV, K_zz_HII_2, c="red")
ax0.fill_between(E/cst.GeV, K_zz_HII_1, K_zz_HII_2, facecolor = "red", alpha = 0.5, label = "HII")

ax0.loglog(E/cst.GeV, K_zz_WIM_1, c="green", ls = "--")
ax0.loglog(E/cst.GeV, K_zz_WIM_2, c="green")
ax0.fill_between(E/cst.GeV, K_zz_WIM_1, K_zz_WIM_2, facecolor = "orange", alpha = 0.5, label = "WIM")

ax0.get_xaxis().set_visible(False)
ax0.legend(loc="upper left")

ax0.set_ylabel("$\\kappa_{zz}$ [cm$^2$s$^{-1}$]")


ax1 = fig.add_subplot(gs[0, 1])

ax1.loglog(E/cst.GeV, K_zz_bc_1, c="black", ls = '--')
ax1.loglog(E/cst.GeV, K_zz_bc_2, c="black")
ax1.fill_between(E/cst.GeV, K_zz_bc_1, K_zz_bc_2, facecolor = "black", alpha = 0., hatch="\\")

ax1.loglog(E/cst.GeV, K_zz_WNM_1, c="green", ls = "--")
ax1.loglog(E/cst.GeV, K_zz_WNM_2, c="green")
ax1.fill_between(E/cst.GeV, K_zz_WNM_1, K_zz_WNM_2, facecolor = "green", alpha = 0.5, label="WNM")

ax1.loglog(E/cst.GeV, K_zz_CNM_1, c="lightblue", ls = "--")
ax1.loglog(E/cst.GeV, K_zz_CNM_2, c="lightblue")
ax1.fill_between(E/cst.GeV, K_zz_CNM_1, K_zz_CNM_2, facecolor = "lightblue", alpha = 0.9, label="CNM")

ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)

ax1.get_xaxis().set_visible(False)
ax1.legend(loc="upper left")

ax2 = fig.add_subplot(gs[1, 0])

ax2.loglog(E/cst.GeV, K_zz_bc_1, c="black", ls = '--')
ax2.loglog(E/cst.GeV, K_zz_bc_2, c="black")
ax2.fill_between(E/cst.GeV, K_zz_bc_1, K_zz_bc_2, facecolor = "black", alpha = 0., hatch="\\")

ax2.loglog(E/cst.GeV, K_zz_DiM_1, c="blue", ls = "--")
ax2.loglog(E/cst.GeV, K_zz_DiM_2, c="blue")
ax2.fill_between(E/cst.GeV, K_zz_DiM_1, K_zz_DiM_2, facecolor = "blue", alpha = 0.5, label = "DiM")

# ax2.get_xaxis().set_visible(False)
ax2.legend(loc="upper left")
ax2.set_ylabel("$\\kappa_{zz}$ [cm$^2$s$^{-1}$]")
ax2.set_xlabel("$E$ [GeV]")

ax3 = fig.add_subplot(gs[1, 1])

ax3.loglog(E/cst.GeV, K_zz_bc_1, c="black", ls = '--')
ax3.loglog(E/cst.GeV, K_zz_bc_2, c="black")
ax3.fill_between(E/cst.GeV, K_zz_bc_1, K_zz_bc_2, facecolor = "black", alpha = 0., hatch="\\")

ax3.loglog(E/cst.GeV, K_zz_DeM_1, c="violet", ls = "--")
ax3.loglog(E/cst.GeV, K_zz_DeM_2, c="violet")
ax3.fill_between(E/cst.GeV, K_zz_DeM_1, K_zz_DeM_2, facecolor = "violet", alpha = 0.5, label = "DeM")

ax3.loglog(E/cst.GeV, K_zz_DeC_1, c="darkblue", ls = "--")
ax3.loglog(E/cst.GeV, K_zz_DeC_2, c="darkblue")
ax3.fill_between(E/cst.GeV, K_zz_DeC_1, K_zz_DeC_2, facecolor = "darkblue", alpha = 0.5, label = "DeC")

ax3.get_yaxis().set_visible(False)
# ax3.get_xaxis().set_visible(False)
ax3.legend(loc="upper left")

ax3.set_xlabel("$E$ [GeV]")


from matplotlib.lines import Line2D

custom_lines = [Line2D([0],[0], color = "white", label = "ISM independant"),
                Line2D([0],[0], color = "white", label = "ISM dependant"),
                Line2D([0],[0], color = "black", ls = '--', label = "$\\delta = 0.5$, $d_{00} = 10^{28}$ cm$^{2}$s$^{-1}$"),
                Line2D([0],[0], color = "black", ls = '--', label = "$q = 3/2$, $I_\\pm = 10^{-4}$"),
                Line2D([0],[0], color = "black", ls = '-', label = "$\\delta = 2/3$, $d_{00} = 10^{29}$ cm$^{2}$s$^{-1}$"),
                Line2D([0],[0], color = "black", ls = '-', label = "$q = 5/3$, $I_\\pm = 10^{-4}$")
                ]


fig.legend(handles=custom_lines, bbox_to_anchor=(0.92, 0.995), ncol = 3)


ymin = 1e27
ymax = 1e33

xmin = 1e-1
xmax = 1e5

ax0.set_ylim(ymin, ymax)
ax1.set_ylim(ymin, ymax)
ax2.set_ylim(ymin, ymax)
ax3.set_ylim(ymin, ymax)

ax0.set_xlim(xmin, xmax)
ax1.set_xlim(xmin, xmax)
ax2.set_xlim(xmin, xmax)
ax3.set_xlim(xmin, xmax)


# plt.plot(E/cst.GeV, K_zz_HII_1, c="red")
# plt.plot(E/cst.GeV, K_zz_WIM_1, c="orange")
# plt.plot(E/cst.GeV, K_zz_WNM_1, c="green")
# plt.plot(E/cst.GeV, K_zz_CNM_1, c="lightblue")
# plt.plot(E/cst.GeV, K_zz_DiM_1, c="blue")
# plt.plot(E/cst.GeV, K_zz_DeM_1, c="violet")
# plt.plot(E/cst.GeV, K_zz_DeC_1, c="black")


# plt.loglog(E/cst.GeV, K_zz_1, c="green")
# plt.loglog(E/cst.GeV, K_zz_2, c="blue")

plt.savefig("diffusion_coeff.pdf")

