#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:20:05 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt
import constants as cst 

# ISM-PHASE FUNCTION 
# This fuction return a table of the 
# different physical properties of the phase 
# T  : Temperature in [K]
# B  : Magnetic field in [G]
# ni : Density of ion [cm^-3]
# X  : Ionization rate 
# mi : Characteristic mass of the ion specie [g]
# mn : Characteristic mass of the neutral specie [g]
# nn : Density of neutral [cm^-3]
# nt : Total density [cm^-3]  
def ism_phase(Temp, Bfiel, nion, ntot, mion, mneutral) : 
    nneutral = ntot - nion
    Xion = nion/ntot
    return dict(T=Temp, B=Bfiel, ni=nion, nn=nneutral, nt=ntot, X=Xion, mi=mion, mn=mneutral)



# Collection of phases 
HII = ism_phase(8000, 10.e-6,  100, 100, cst.mHII, 0.93*cst.mHI+0.07*cst.mHeI)
WIM = ism_phase(8000, 5.00001e-6,  0.315, 0.35, cst.mHII, 0.93*cst.mHI+0.07*cst.mHeI)
WNM = ism_phase(8000, 5.00001e-6,   7e-3, 0.35, cst.mHII, 0.93*cst.mHI+0.07*cst.mHeI)
CNM = ism_phase(  50, 6.00001e-6, 2.3e-2, 30.0, cst.mCII, 0.93*cst.mHI+0.07*cst.mHeI)
DiM = ism_phase(  50, 6.00001e-6, 3.0e-2,  300, cst.mCII, 0.93*(0.5*cst.mHI + 0.5*cst.mH2) + 0.07*cst.mHeI)

