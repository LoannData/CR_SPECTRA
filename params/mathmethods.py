#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:04:53 2018

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 





# Methode analytique de cardano (x**3 + a x**2 + b x + c = 0)
# Qui renvoie une solution réelle uniquement
def cardano3(a, b, c) : 
    p = b - a**2/3.
    q = 2*a**3/27. - a*b/3. + c
    R = q/2.
    Q = p/3.
    D = Q**3 + R**2
    x = 0.
    if (D >= 0.) : 
        if (-R + np.sqrt(D) >= 0.):
            S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R + np.sqrt(D) < 0.):
            S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R - np.sqrt(D) >= 0.): 
            S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
        if (-R -np.sqrt(D) < 0.) : 
            S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
        x = -(1/3.)*a + (S1 + S2) #Solution réelle 
    if (D < 0.) : 
        D = - D
        if (-R + np.sqrt(D) >= 0.):
            S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R + np.sqrt(D) < 0.):
            S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R - np.sqrt(D) >= 0.): 
            S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
        if (-R -np.sqrt(D) < 0.) : 
            S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
        x = -(1/3.)*a + (S1 + S2) #Solution réelle 
    return x 


def histogram(data, xi, xf, nbin, scale, normalization) : 
    if (scale == 'linear') : 
#        if (not(xi)) : xi = min(data)
#        if (not(xf)) : xf = max(data)
        delta_x = (float(xf) - float(xi))/nbin
        xnew = np.linspace(xi+delta_x/2., xf-delta_x/2., nbin)
        distribution = np.zeros(len(xnew))
        for ii in range(len(xnew)) : 
            for jj in range(len(data)) : 
                if (data[jj] < xnew[ii]+delta_x/2. and data[jj] >= xnew[ii]-delta_x/2.) : 
                    
                    distribution[ii] += 1
        if (normalization) : 
#            distribution = distribution/len(data)*normalization
            distribution = distribution*normalization
        return xnew, distribution



