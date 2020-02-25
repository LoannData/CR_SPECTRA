#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:21:49 2018

@author: lbrahimi
"""

import numpy as np

# Constants 
yr   = 3.154e+7   # 1 yr in [s]
kyr  = 1e3*yr     # 1 kyr in [s]
pc   = 3.086e18   # 1 pc in [cm]
kpc  = 1.e3*pc    # 1 kpc in [cm]
GeV  = 0.00160218 # 1 GeV in [erg] 
TeV  = 1e3*GeV    # 1 TeV in [erg]
eV   = GeV*1e-9   # 1 eV in [erg]
MeV  = 1e-3*GeV
me   = 9.10938e-28   # Electron mass [g]
mp   = 1.6726219e-24 # Proton mass [g]
mn   = 1.6749286e-24 # Neutron mass [g]
mHI  = mp # HI mass [g]
mHII  = mp # HII mass [g]
mHeI = 2*mp + 2*mn # HeI mass [g]
mCII = 6*mp + 6*mn # CII mass [g]
mH2  = 2*mp        # H2 mass [g]
e = 4.80326e-10 # e charge in [statC]
c = 29979245800. # light celerity in [cm/s]
kbolz = 1.3807e-16 # Boltzmann constant in CGS
kms = 1e5 # 1 km/s in cm/s 