#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:28:36 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 


# import constants as cst

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
    if (name == "symlog") : 
        S = np.empty(Ns)
        S_right = np.empty(int(Ns/2.))
        S_left  = np.empty(int(Ns/2.))
        S_middle = (Smax - Smin)/2. 
        loc_S = np.logspace(np.log10(Smin), np.log10(S_middle), int(Ns/2.))
        
        
        for si in range(int(Ns/2.)) : 
            S_right[si] = S_middle + loc_S[si] - Smin 
        for si in range(int(Ns/2.)) : 
            S_left[si] = S_middle - (loc_S[si] - Smin)
            
        for si in range(int(Ns/2.)) : 
            S[si] = S_left[-1-si]
        for si in range(int(Ns/2.), Ns) : 
            S[si] = S_right[si - int(Ns/2.)]
        return S

    
# X = grid(100.*cst.pc, 1000*cst.pc, 2**8, "symlog")
        
        
        
        
