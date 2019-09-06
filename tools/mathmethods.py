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




# Methode permettant d'adoucir la transition entre deux phases suivant une 
# largeur bien définie 
def g1(x, xt, l) : 
    return 0.5*(1 + np.tanh((xt-x)/l))
#    return 0.

def g2(x, xt, l) : 
    return 0.5*(1 + np.tanh(-(xt-x)/l))
#    return 0.

def f(x, xt, l, v1, v2) : 
    return g1(x, xt, l)*v1 + g2(x, xt, l)*v2

def SmoothPhaseTransition(X, phases, smooth_width) : 
    T  = np.zeros(len(X))
    B  = np.zeros(len(X))
    ni = np.zeros(len(X))
    nn = np.zeros(len(X))
    nt = np.zeros(len(X))
    Xi  = np.zeros(len(X))
    mi = np.zeros(len(X))
    mn = np.zeros(len(X))
    if (len(phases) == 1) : 
        for xi in range(len(X)) : 
            for pi in range(len(phases)) : 
                if (X[xi] >= phases[pi][1].get("Xmin") and X[xi] <= phases[pi][1].get("Xmax")) : 
                    T[xi]  = phases[pi][0].get("T")
                    B[xi]  = phases[pi][0].get("B")
                    ni[xi] = phases[pi][0].get("ni")
                    nn[xi] = phases[pi][0].get("nn")
                    nt[xi] = phases[pi][0].get("nt")
                    Xi[xi] = phases[pi][0].get("X")
                    mi[xi] = phases[pi][0].get("mi")
                    mn[xi] = phases[pi][0].get("mn")
    if (len(phases) > 1) : 
        for xi in range(len(X)) :
            for pos in range(len(phases)) :        
                if (X[xi] >= phases[pos][1].get('Xmin') and X[xi] < phases[pos][1].get('Xmin') + 0.5*(phases[pos][1].get('Xmax') - phases[pos][1].get('Xmin'))) : 
                    xt = phases[pos][1].get('Xmin')
                    l  = smooth_width
                    if (pos == 0) :
                        v1 = phases[pos][0]
                        v2 = phases[pos][0]
                    else : 
                        v1 = phases[pos-1][0]
                        v2 = phases[pos][0]
                    T[xi]  = f(X[xi], xt, l, v1.get('T'), v2.get('T'))
                    B[xi]  = f(X[xi], xt, l, v1.get('B'), v2.get('B'))
                    ni[xi]  = f(X[xi], xt, l, v1.get('ni'), v2.get('ni'))
                    nn[xi]  = f(X[xi], xt, l, v1.get('nn'), v2.get('nn'))
                    nt[xi] = f(X[xi], xt, l, v1.get('nt'), v2.get('nt'))
                    Xi[xi]  = f(X[xi], xt, l, v1.get('X'), v2.get('X'))
                    mi[xi]  = f(X[xi], xt, l, v1.get('mi'), v2.get('mi'))
                    mn[xi]  = f(X[xi], xt, l, v1.get('mn'), v2.get('mn'))
                if (X[xi] >= phases[pos][1].get('Xmin') + 0.5*(phases[pos][1].get('Xmax') - phases[pos][1].get('Xmin')) and X[xi] < phases[pos][1].get('Xmax')) :
                    xt = phases[pos][1].get('Xmax')
                    l  = smooth_width
                    if (pos == len(phases)-1) :
                        v1 = phases[pos][0]
                        v2 = phases[pos][0]
                    else : 
                        v1 = phases[pos][0]
                        v2 = phases[pos+1][0]
                    T[xi]  = f(X[xi], xt, l, v1.get('T'), v2.get('T'))
                    B[xi]  = f(X[xi], xt, l, v1.get('B'), v2.get('B'))
                    ni[xi]  = f(X[xi], xt, l, v1.get('ni'), v2.get('ni'))
                    nn[xi]  = f(X[xi], xt, l, v1.get('nn'), v2.get('nn'))
                    nt[xi] = f(X[xi], xt, l, v1.get('nt'), v2.get('nt'))
                    Xi[xi]  = f(X[xi], xt, l, v1.get('X'), v2.get('X'))
                    mi[xi]  = f(X[xi], xt, l, v1.get('mi'), v2.get('mi'))
                    mn[xi]  = f(X[xi], xt, l, v1.get('mn'), v2.get('mn'))
    return T, B, ni, nn, nt, Xi, mi, mn


