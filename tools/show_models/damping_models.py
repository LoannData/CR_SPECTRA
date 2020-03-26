#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:19:40 2020

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import sys 

sys.path.append("../")
import constants as cst
import phases_collection as ism
import damping as dp
import d1_grid_generator as grid 

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    
    




def damprate_to_damptime(gamma) : 
    return gamma**(-1)/cst.yr 
        
        
        












# We create the energy grid 
NE        = 10  # 2**NE is the E-size of the grid
Emin      = 0.99*cst.GeV
Emax      = 100.01*cst.TeV
egridtype = "logspace" # Type of grid - # logspace type recomended (only option for the moment)
E = grid.grid(Emin, Emax, 2**NE, egridtype)  




phases = [ism.HII, ism.WIM, ism.WNM, ism.CNM, ism.DiM, ism.DeM, ism.DeC]

xlim   = [Emin/cst.GeV, Emax/cst.GeV]
xlims  = [xlim, xlim, xlim, xlim, xlim, xlim, xlim]

ylim  = [1e-14, 1e-3]
ylims = [ylim, ylim, ylim, ylim, ylim, ylim, ylim]

Name = ["HII", "WIM", "WNM", "CNM", "DiM", "DeM", "DeC"]

size_x = 4
size_y = 3
sub_x  = 2
sub_y  = 4
fig = plt.figure(figsize=(size_x*sub_x,size_y*sub_y))

gs = gridspec.GridSpec(ncols= sub_x, nrows = sub_y, figure = fig )
gs.update(wspace=0.05, hspace=0.05) # set the spacing between axes.


pos_1 = [0, 0, 1, 1, 2, 2, 3, 3]
pos_2 = [0, 1, 0, 1, 0, 1, 0, 1]

for pi in range(len(phases)) : 
    # Alfven waves propagation
    wR_Alfven = np.zeros(len(E))
    wI_Alfven = np.zeros(len(E))
    
    wR_Alfven_o1 = np.zeros(len(E))
    wI_Alfven_o1 = np.zeros(len(E))
    
    wR_Alfven_o2 = np.zeros(len(E))
    wI_Alfven_o2 = np.zeros(len(E))
    
    Ep = np.NaN
    Em = np.NaN
    
    
    Gamma_lz     = np.zeros(len(E))
    
    Gamma_nlld_inf   = np.zeros(len(E))
    Gamma_nlld_sup   = np.zeros(len(E))

    
    
    for e in range(len(E)) : 
        
        in_damping = dp.IonNeutral_Damping(E[e], phases[pi], nu_n = 0, theta = 0)
        wR_Alfven[e] = in_damping.get("wr")
        wI_Alfven[e] = in_damping.get("wi")
        
        if (pi == 0) : 
            wI_Alfven[e] = np.NaN
        
        
        in_damping = dp.IN_damping_approx_1(E[e], phases[pi], nu_n = 0, theta = 0)
        wR_Alfven_o1[e] = in_damping.get("wr")
        wI_Alfven_o1[e] = in_damping.get("wi")
        
        in_damping = dp.IN_damping_approx_2(E[e], phases[pi], theta = 0)
        wR_Alfven_o2[e] = in_damping.get("wr")
        wI_Alfven_o2[e] = in_damping.get("wi")
        Ep = in_damping.get("Ep")
        Em = in_damping.get("Em")
        
        if (pi > 1) : 
            lz_damping = dp.damping_lazarian_nopos(E[e], phases[pi])
            Gamma_lz[e] = lz_damping[0]
            lz_min      = lz_damping[1]
        
        Iinf = 1e-4 
        Isup = 1e-1
        Gamma_nlld_inf[e] = dp.non_linear_landau_damping(phases[pi].get("T"), Iinf, Iinf, 
                                                         phases[pi].get("mi"), cst.e, phases[pi].get("B"), E[e])
        Gamma_nlld_sup[e] = dp.non_linear_landau_damping(phases[pi].get("T"), Isup, Isup, 
                                                         phases[pi].get("mi"), cst.e, phases[pi].get("B"), E[e])
        
        
        
    
    # plt.figure(figsize = (12, 8))
    ax = fig.add_subplot(gs[pos_1[pi], pos_2[pi]])
    
    ax.loglog(E/cst.GeV, wR_Alfven, c="blue")
    ax.loglog(E/cst.GeV, -wI_Alfven, c="black")
    
    ax.loglog(E/cst.GeV, wR_Alfven_o1, c="blue", ls='--')
    ax.loglog(E/cst.GeV, wI_Alfven_o1, c="black", ls='--')
    
    ax.loglog(E/cst.GeV, wR_Alfven_o2, c="blue", ls=':', lw = 5)
    ax.loglog(E/cst.GeV, wI_Alfven_o2, c="black", ls=':', lw = 5)
    
    # plt.axvline(Ep/cst.GeV, c="black")
    # plt.axvline(Em/cst.GeV, c="black")
    ax.axvspan(Em/cst.GeV, Ep/cst.GeV, alpha=0.5, color='grey')
    
    if (pi > 1) : 
        ax.loglog(E/cst.GeV, -Gamma_lz, c="red", ls = '-')
        ax.axvline(lz_min/cst.GeV, c="red", lw = 5)
    
    # plt.loglog(E/cst.GeV, -Gamma_nlld_inf, c="orange")
    # plt.loglog(E/cst.GeV, -Gamma_nlld_sup, c="orange")
    ax.fill_between(E/cst.GeV, -Gamma_nlld_inf, -Gamma_nlld_sup, facecolor = "yellow", alpha = 0.5)
    
    ax.set_xlim(xlims[pi][0], xlims[pi][1])
    ax.set_ylim(ylims[pi][0], ylims[pi][1])
    
    if (pos_1[pi] != 3 and not(pos_1[pi] == 2 and pos_2[pi] == 1)) : 
        ax.get_xaxis().set_visible(False)
    else : 
        ax.set_xlabel("E [GeV]")
    if (pos_2[pi] == 1) : 
        ax.get_yaxis().set_visible(False)
    else : 
        ax.set_ylabel("$-\\Gamma, \\omega$ [s$^{-1}]$")
    
    ax.set_title(Name[pi], x = 0.9, y = 0.8)
    
    if (pos_1[pi] == 3 and pos_2[pi] == 0) : 
        ax.plot([],[], c = "blue", label = "$\\omega_R^A$ Exact solution")
        ax.plot([],[], c = "black", label = "$-\\Gamma_I^A$ Exact solution")
        ax.plot([],[], c = "blue", ls = '--', label = "$\\omega_R^A$ where $-\\Gamma_I^A \\ll \\omega_R^A$")
        ax.plot([],[], c = "black", ls = '--', label = "$-\\Gamma_I^A$ where $-\\Gamma_I^A \\ll \\omega_R^A$")
        ax.plot([],[], c = "red", ls = "-", label = "L16 damping model")
        ax.fill_between([],[],[], facecolor="yellow", alpha = 0.5, label="NLLD for $10^{-4} < \\delta B^2/B_0^2 < 10^{-1}$")
        ax.legend(loc = 1, bbox_to_anchor=(2.,0.8)) 

# plt.tight_layout()
plt.savefig("damping_terms.pdf", bbox_inches='tight',pad_inches=-0.)