#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:17:05 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import sys 
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from scipy import interpolate


# We refers to the tools setup folder 
import freader as fr 
import fwritter as fw
import constants as cst
import d1_grid_generator as grid 
import phases_collection as ism
import damping as dp

nx        = int(fr.search("gparams.dat","NX"))
ne        = int(fr.search("gparams.dat","NE")) 



# Size of the problem 
#nx        = 12
NX        = 2**nx # X position array 
Xmin      = 0.*cst.pc # Lower boundary position of the box in [cm]
Xmax      = 500.*cst.pc # Upper boundary position of the box in [cm]
xgridtype = "cartesian" # Type of grid 

#ne        = 0
NE        = 2**ne # Energy array 
Emin      = 9.99*cst.GeV # Minimum energy of the distribution [erg]
Emax      = 10.01*cst.TeV # Maximum energy of the distribution [erg]
egridtype = "logspace" # Type of grid - # logspace type recomended 


# Phase space 
X = grid.grid(Xmin, Xmax, NX, xgridtype)
E = grid.grid(Emin, Emax, NE, egridtype)  

# ISM composition 
phases = []
phases.append([ism.WNM, dict(Xmin=0.*cst.pc, Xmax=500.*cst.pc)])


T  = np.zeros(len(X))
B  = np.zeros(len(X))
ni = np.zeros(len(X))
nn = np.zeros(len(X))
nt = np.zeros(len(X))
Xi  = np.zeros(len(X))
mi = np.zeros(len(X))
mn = np.zeros(len(X))

smooth_width = 3.*cst.pc # Smooth width 


def g1(x, xt, l) : 
    return 0.5*(1 + np.tanh((xt-x)/l))
#    return 0.

def g2(x, xt, l) : 
    return 0.5*(1 + np.tanh(-(xt-x)/l))
#    return 0.

def f(x, xt, l, v1, v2) : 
    return g1(x, xt, l)*v1 + g2(x, xt, l)*v2

for xi in range(len(X)) : 
    T[xi]  = phases[0][0].get("T")
    B[xi]  = phases[0][0].get("B")
    ni[xi] = phases[0][0].get("ni")
    nn[xi] = phases[0][0].get("nn")
    nt[xi] = phases[0][0].get("nt")
    Xi[xi] = phases[0][0].get("X")
    mi[xi] = phases[0][0].get("mi")
    mn[xi] = phases[0][0].get("mn")

Xcr_min = 0.*cst.pc
Xcr_max = 20.*cst.pc#24.72*cst.pc

smooth_cr_width = (Xcr_max-Xcr_min)/100.
alfven_width = smooth_cr_width#3*cst.pc
        
# ISM secondary variables 
ism_values = dict(T=T, B=B, ni=ni, nn=nn, nt=nt, X=Xi, mi=mi, mn=mn)

# Calculation of the Alfven velocity and the ion-neutral damping 
VA       = np.zeros((len(E), len(X)))
gamma_in  = np.zeros((len(E), len(X)))
gamma_lazarian = np.zeros((len(E), len(X)))
gamma_tot = np.zeros((len(E), len(X)))
Pcr       = np.zeros((len(E), len(X)))
for e in range(len(E)) : 
    print (e,"/",len(E))
    for xi in range(len(X)) : 
        in_damping = dp.indamping_alfven(xi , E[e], ism_values) 
        VA[e][xi] = in_damping.get('VA')*1e0#f(X[xi], Xcr_max, alfven_width, 0., in_damping.get('VA'))
        gamma_in[e][xi] = in_damping.get('wi')
        gamma_lazarian[e][xi] = -dp.damping_lazarian(xi, E[e], ism_values)
        gamma_tot[e][xi] = gamma_in[e][xi] + gamma_lazarian[e][xi]
        Pcr[e][xi] = 0.0001*cst.eV

# Initial diffusion coefficients -> Initial waves I
d00     = np.zeros(len(X)) # Spatial component of the background diffusion coefficient of CRs 
D       = np.zeros((len(E), len(X))) # Initial diffusion coefficient of CRs 
Db      = np.zeros((len(E), len(X))) # Bohm diffusion coefficient of CRs 
I       = np.zeros((len(E), len(X))) # Initial resonnant waves energy density 


for xi in range(len(X)) : 
    d00[xi] = 1e28


for e in range(len(E)) : 
    print (e,"/",len(E))
    for xi in range(len(X)) : 
        Db[e][xi] = (4*cst.c)/(3*np.pi)*(E[e]/(cst.e*B[xi]))
        D[e][xi]  = d00[xi]*(E[e]/(10.*cst.GeV))**0.5
        I[e][xi] = Db[e][xi]/D[e][xi]

fw.write1D( "Alfven.dat", nx=nx, ne=ne, variable=VA.T, path="../data_ini/") 
fw.write1D("DCRpara.dat", nx=nx, ne=ne,  variable=D.T, path="../data_ini/")
fw.write1D(  "DBohm.dat", nx=nx, ne=ne, variable=Db.T, path="../data_ini/")
fw.write1D(     "Ip.dat", nx=nx, ne=ne,  variable=I.T, path="../data_ini/")   
fw.write1D(    "Pcr.dat", nx=nx, ne=ne,variable=Pcr.T, path="../data_ini/")    
fw.write1D("damping.dat", nx=nx, ne=ne, variable=gamma_tot.T, path="../data_ini/")  

fw.write1Daxis("X.dat", variable=X, nx=nx, path="../data_ini/")
fw.write1Daxis("E.dat", variable=E, nx=ne, path="../data_ini/")

fw.write1Daxis("B.dat", variable=B, nx=nx, path="../data_ini/")


# We write some important informations about SNRs and CRs 
Esn = 1e51 # [erg]
Mej = 1    # ejecta mass in sum mass units 
xi_n = 1   # for solar abundances 
phi_c = 1 # Actual thermal cond. / the Sptitzer (1962) value 
beta  = 2. 


variables = {"ni"      : ni[0],
             "X"       : Xi[0],
             "mn"      : mn[0],
             "T"       : T[0], 
             "Esn"     : Esn, 
             "Mej"     : Mej,
             "xi_n"    : xi_n,
             "phi_c"   : phi_c, 
             "beta"    : beta}

fw.fileWrite("specparams", variables = variables, path='./', ext='.dat') 
















  