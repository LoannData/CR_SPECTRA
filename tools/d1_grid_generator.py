#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:28:36 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import scipy.special as sp

import constants as cst 

def grid(Smin, Smax, Ns, name, s_center = None, width = None, smooth = None, dXmin = 0.01*cst.pc) : 
    if (name == "cartesian") : 
        S = np.empty(Ns)
#        S = np.linspace(Smin, Smax, Ns+2)
        for si in range(Ns) : 
            S[si] = Smin + (Smax - Smin)/(Ns)*si
        return S
    if (name == "logspace") : 
        S = np.empty(Ns)
        S = np.logspace(np.log10(Smin), np.log10(Smax), Ns)
        return S
    if (name == "gaussian") : # Not really working 
        mu = s_center
        sig = 60.*cst.pc 
        r = 0.001 
        alpha = 1 - r
        
        C1 = (Smax/cst.pc - Smin/cst.pc)
        C2 = (Smax/cst.pc - Smin/cst.pc)
        C3 = - 1/np.sqrt(2)*alpha*sig/cst.pc*(sp.erf((mu - Smin)/sig) - sp.erf((mu - Smax)/sig))
        C = C1/(C2 + C3)/(Ns)*(Smax - Smin) 
        

        S = grid(Smin, Smax, Ns, "cartesian") 
        
        dX = np.empty(Ns)
        for xi in range(Ns) : 
            dX[xi] = C*(1 - alpha*np.exp(-(S[xi] - mu)**2/sig**2))
        
        X = np.empty(Ns)
        X[0] = Smin 
        for xi in range(1, Ns) : 
            X[xi] = X[xi-1] + dX[xi] 
            
            
        # plt.plot(dX/cst.pc)  
        # plt.plot(X/cst.pc)
        
        print ("Number of bins : ",Ns)
        print ("Normal step = ",(Smax - Smin)/Ns/cst.pc," pc")
        print ("Inverted gaussian parameters :")
        print ("mean = ",mu/cst.pc," pc, sigma = ",sig/cst.pc," pc")
        print ("Ratio small/large bins : ",r)
        print ("min Step = ",min(dX)/cst.pc," pc")
        print ("max Step = ",max(dX)/cst.pc," pc")
        print ("Xmin = ", X[0]/cst.pc," pc, Xmax = ",X[-1]/cst.pc," pc")
        return X
    if (name == "erf") : # Not working yet 
        mu = width/2.
        sig = smooth 
        dXmin = dXmin

        S = grid(Smin, Smax, Ns, "cartesian") 
        
        dX = np.empty(Ns)
        for xi in range(Ns) : 
            dX[xi] = dXmin + 0.5*(Smax - Smin + width)/Ns*(2 + (sp.erf((S[xi]-s_center-mu)/sig) - sp.erf((S[xi]-s_center+mu)/sig)))
        
        X = np.empty(Ns)
        X[0] = Smin 
        for xi in range(1, Ns) : 
            X[xi] = X[xi-1] + dX[xi] 
        
        print ("Number of bins : ",Ns)
        print ("Normal step = ",(Smax - Smin)/Ns/cst.pc," pc")
        print ("Erf parameters : ")
        print ("Center position : ",s_center/cst.pc," pc, width : ",width/cst.pc," pc, smooth : ",smooth/cst.pc," pc")
        print ("min DX = ",dXmin/cst.pc," pc, max DX = ",max(dX)/cst.pc," pc")
        
        return X
        

        
        




# X = grid(0.*cst.pc, 2000.*cst.pc, 2**9, "erf", 
#          s_center = 1000.*cst.pc, width = 200.*cst.pc, smooth = 20.*cst.pc) 