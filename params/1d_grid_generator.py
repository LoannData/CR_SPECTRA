#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:28:36 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 

def grid(Smin, Smax, Ns, name) : 
    if (name == "cartesian") : 
        S = np.empty(Ns)
#        S = np.linspace(Smin, Smax, Ns+2)
        for si in range(Ns) : 
            S[si] = Smin + (Smax - Smin)/(Ns)*(si+1)
        return S
    if (name == "logspace") : 
        S = np.empty(Ns)
        S = np.logspace(np.log10(Smin), np.log10(Smax), Ns)
        return S