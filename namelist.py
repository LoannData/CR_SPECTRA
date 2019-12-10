#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 08:10:35 2019

@author: lbrahimi
"""
#import numpy as np
#import matplotlib.pyplot as plt 
#from matplotlib.patches import Rectangle
#import matplotlib.gridspec as gridspec
#from scipy import interpolate


# We refers to the tools setup folder
import sys
sys.path.append('./tools/')
 
#import freader as fr 
#import fwritter as fw
import constants as cst
import d1_grid_generator as grid 
import phases_collection as ism
import damping as dmp
#import damping as dp
import mathmethods as mh


# Function to get a mean VA in a given phase 
def getVA(E, phase) : 
    VA = []
    for ii in range(2**NE) : 
        VA.append(dmp.indamping_alfven_nopos(E[ii], phase).get("VA"))    
    return VA

###############################################################################
#      GRID PARAMETERS                                                        #
###############################################################################
NX = 13 # 2**NX is the X-size of the grid
NE = 6  # 2**NE is the E-size of the grid 

Xmin = 0.*cst.pc
Xmax = 2000.*cst.pc
xgridtype = "cartesian"

Emin = 0.99*cst.GeV
Emax = 100.01*cst.TeV
egridtype = "logspace" # Type of grid - # logspace type recomended 

box_center = 1000.*cst.pc 

# Phase space 
X = grid.grid(Xmin, Xmax, 2**NX, xgridtype)
E = grid.grid(Emin, Emax, 2**NE, egridtype)  

###############################################################################
#       OTHER TERMS                                                           #
###############################################################################
in_damping       = True # Ion neutral damping of waves
lz_damping       = True # Lazarian damping of waves
Pcr_1GeV         = 1*cst.eV # [cm^-3] CR Pressure at 1 GeV 


###############################################################################
#        ISM STRUCTURE                                                        #
###############################################################################
phases  = [] # Phases list
# Append phases in the order of the setup you want to create
#phases.append([ism.WNM, dict(Xmin=0.*cst.pc, Xmax=2000.*cst.pc), getVA(E, ism.WNM)]) 
#phases.append([ism.CNM, dict(Xmin=100.*cst.pc, Xmax=200.*cst.pc), getVA(E, ism.CNM)])
#phases.append([ism.DiM, dict(Xmin=510.*cst.pc, Xmax=1000.*cst.pc), getVA(E, ism.DiM)])
#phases.append([ism.CNM, dict(Xmin=230.*cst.pc, Xmax=330.*cst.pc), getVA(E, ism.CNM)])
#phases.append([ism.WNM, dict(Xmin=330.*cst.pc, Xmax=1000.*cst.pc), getVA(E, ism.WNM)])

phases.append([ism.WNM, dict(Xmin=0.*cst.pc, Xmax=300.*cst.pc), getVA(E, ism.WNM)]) 
phases.append([ism.CNM, dict(Xmin=300.*cst.pc, Xmax=500.*cst.pc), getVA(E, ism.CNM)])
phases.append([ism.DiM, dict(Xmin=500.*cst.pc, Xmax=600.*cst.pc), getVA(E, ism.DiM)]) 
phases.append([ism.CNM, dict(Xmin=600.*cst.pc, Xmax=800.*cst.pc), getVA(E, ism.CNM)])
phases.append([ism.WNM, dict(Xmin=800.*cst.pc, Xmax=1200.*cst.pc), getVA(E, ism.WNM)]) 
phases.append([ism.CNM, dict(Xmin=1200.*cst.pc, Xmax=1400.*cst.pc), getVA(E, ism.CNM)]) 
phases.append([ism.DiM, dict(Xmin=1400.*cst.pc, Xmax=1500.*cst.pc), getVA(E, ism.DiM)]) 
phases.append([ism.CNM, dict(Xmin=1500.*cst.pc, Xmax=1700.*cst.pc), getVA(E, ism.CNM)]) 
phases.append([ism.WNM, dict(Xmin=1700.*cst.pc, Xmax=2000.*cst.pc), getVA(E, ism.WNM)]) 

smooth_width_transition = 10.*cst.pc # Smooth width transition between two phases

# We calculate the smoothed variables
T, B, ni, nn, nt, Xi, mi, mn, va = mh.SmoothPhaseTransition(X, E, phases, smooth_width_transition)
# ISM secondary variables 
ism_values = dict(T=T, B=B, ni=ni, nn=nn, nt=nt, X=Xi, mi=mi, mn=mn, VA=va)


