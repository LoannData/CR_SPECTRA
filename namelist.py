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
#import damping as dp
import mathmethods as mh



###############################################################################
#      GRID PARAMETERS                                                        #
###############################################################################
NX = 12 # 2**NX is the X-size of the grid
NE = 6  # 2**NE is the E-size of the grid 

Xmin = 0.*cst.pc
Xmax = 2000.*cst.pc
xgridtype = "cartesian"

Emin = 9.99*cst.GeV
Emax = 10.01*cst.TeV
egridtype = "logspace" # Type of grid - # logspace type recomended 

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
phases = [] # Phases list
# Append phases in the order of the setup you want to create
phases.append([ism.WNM, dict(Xmin=0.*cst.pc, Xmax=100.*cst.pc)]) 
phases.append([ism.CNM, dict(Xmin=100.*cst.pc, Xmax=2000.*cst.pc)])
smooth_width_transition = 10.*cst.pc # Smooth width transition between two phases
# We calculate the smoothed variables
T, B, ni, nn, nt, Xi, mi, mn = mh.SmoothPhaseTransition(X, phases, smooth_width_transition)
# ISM secondary variables 
ism_values = dict(T=T, B=B, ni=ni, nn=nn, nt=nt, X=Xi, mi=mi, mn=mn)


