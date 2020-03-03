#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:04:53 2018

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import math


# Main Function takes in the coefficient of the Cubic Polynomial
# as parameters and it returns the roots in form of numpy array.
# Polynomial Structure -> ax^3 + bx^2 + cx + d = 0

def Cubic3(a, b, c, d):

    if (a == 0 and b == 0):                     # Case for handling Liner Equation
        return np.array([(-d * 1.0) / c])                 # Returning linear root as numpy array.

    elif (a == 0):                              # Case for handling Quadratic Equations

        D = c * c - 4.0 * b * d                       # Helper Temporary Variable
        if D >= 0:
            D = math.sqrt(D)
            x1 = (-c + D) / (2.0 * b)
            x2 = (-c - D) / (2.0 * b)
        else:
            D = math.sqrt(-D)
            x1 = (-c + D * 1j) / (2.0 * b)
            x2 = (-c - D * 1j) / (2.0 * b)
            
        return np.array([x1, x2])               # Returning Quadratic Roots as numpy array.

    f = findF(a, b, c)                          # Helper Temporary Variable
    g = findG(a, b, c, d)                       # Helper Temporary Variable
    h = findH(g, f)                             # Helper Temporary Variable

    if f == 0 and g == 0 and h == 0:            # All 3 Roots are Real and Equal
        if (d / a) >= 0:
            x = (d / (1.0 * a)) ** (1 / 3.0) * -1
        else:
            x = (-d / (1.0 * a)) ** (1 / 3.0)
        return np.array([x, x, x])              # Returning Equal Roots as numpy array.

    elif h <= 0:                                # All 3 roots are Real

        i = math.sqrt(((g ** 2.0) / 4.0) - h)   # Helper Temporary Variable
        j = i ** (1 / 3.0)                      # Helper Temporary Variable
        k = math.acos(-(g / (2 * i)))           # Helper Temporary Variable
        L = j * -1                              # Helper Temporary Variable
        M = math.cos(k / 3.0)                   # Helper Temporary Variable
        N = math.sqrt(3) * math.sin(k / 3.0)    # Helper Temporary Variable
        P = (b / (3.0 * a)) * -1                # Helper Temporary Variable

        x1 = 2 * j * math.cos(k / 3.0) - (b / (3.0 * a))
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P

        return np.array([x1, x2, x3])           # Returning Real Roots as numpy array.

    elif h > 0:                                 # One Real Root and two Complex Roots
        R = -(g / 2.0) + math.sqrt(h)           # Helper Temporary Variable
        if R >= 0:
            S = R ** (1 / 3.0)                  # Helper Temporary Variable
        else:
            S = (-R) ** (1 / 3.0) * -1          # Helper Temporary Variable
        T = -(g / 2.0) - math.sqrt(h)
        if T >= 0:
            U = (T ** (1 / 3.0))                # Helper Temporary Variable
        else:
            U = ((-T) ** (1 / 3.0)) * -1        # Helper Temporary Variable

        x1 = (S + U) - (b / (3.0 * a))
        x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * math.sqrt(3) * 0.5j
        x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * math.sqrt(3) * 0.5j

        return np.array([x1, x2, x3])           # Returning One Real Root and two Complex Roots as numpy array.


# Helper function to return float value of f.
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0


# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0


# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)





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

def SmoothPhaseTransition(X, E, phases, smooth_width) : 
    T  = np.zeros(len(X))
    B  = np.zeros(len(X))
    ni = np.zeros(len(X))
    nn = np.zeros(len(X))
    nt = np.zeros(len(X))
    Xi  = np.zeros(len(X))
    mi = np.zeros(len(X))
    mn = np.zeros(len(X))
    VA = np.zeros((len(E), len(X)))
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
                    for ei in range(len(E)) : 
                        VA[ei][xi] = phases[pi][2][ei]
    if (len(phases) > 1) : 
        for xi in range(len(X)) :
            for pos in range(len(phases)) :        
                if (X[xi] >= phases[pos][1].get('Xmin') and X[xi] < phases[pos][1].get('Xmin') + 0.5*(phases[pos][1].get('Xmax') - phases[pos][1].get('Xmin'))) : 
                    xt = phases[pos][1].get('Xmin')
                    l  = smooth_width
                    if (pos == 0) :
                        v1 = phases[pos][0]
                        v2 = phases[pos][0]
                        va1 = phases[pos][2]
                        va2 = phases[pos][2]
                    else : 
                        v1 = phases[pos-1][0]
                        v2 = phases[pos][0]
                        va1 = phases[pos-1][2]
                        va2 = phases[pos][2]
                    T[xi]  = f(X[xi], xt, l, v1.get('T'), v2.get('T'))
                    B[xi]  = f(X[xi], xt, l, v1.get('B'), v2.get('B'))
                    ni[xi]  = f(X[xi], xt, l, v1.get('ni'), v2.get('ni'))
                    nn[xi]  = f(X[xi], xt, l, v1.get('nn'), v2.get('nn'))
                    nt[xi] = f(X[xi], xt, l, v1.get('nt'), v2.get('nt'))
                    Xi[xi]  = f(X[xi], xt, l, v1.get('X'), v2.get('X'))
                    mi[xi]  = f(X[xi], xt, l, v1.get('mi'), v2.get('mi'))
                    mn[xi]  = f(X[xi], xt, l, v1.get('mn'), v2.get('mn'))
                    for ei in range(len(E)) : 
                        VA[ei][xi] = f(X[xi], xt, l, va1[ei], va2[ei])
                if (X[xi] >= phases[pos][1].get('Xmin') + 0.5*(phases[pos][1].get('Xmax') - phases[pos][1].get('Xmin')) and X[xi] < phases[pos][1].get('Xmax')) :
                    xt = phases[pos][1].get('Xmax')
                    l  = smooth_width
                    if (pos == len(phases)-1) :
                        v1 = phases[pos][0]
                        v2 = phases[pos][0]
                        va1 = phases[pos][2]
                        va2 = phases[pos][2]
                    else : 
                        v1 = phases[pos][0]
                        v2 = phases[pos+1][0]
                        va1 = phases[pos][2]
                        va2 = phases[pos+1][2]
                    T[xi]  = f(X[xi], xt, l, v1.get('T'), v2.get('T'))
                    B[xi]  = f(X[xi], xt, l, v1.get('B'), v2.get('B'))
                    ni[xi]  = f(X[xi], xt, l, v1.get('ni'), v2.get('ni'))
                    nn[xi]  = f(X[xi], xt, l, v1.get('nn'), v2.get('nn'))
                    nt[xi] = f(X[xi], xt, l, v1.get('nt'), v2.get('nt'))
                    Xi[xi]  = f(X[xi], xt, l, v1.get('X'), v2.get('X'))
                    mi[xi]  = f(X[xi], xt, l, v1.get('mi'), v2.get('mi'))
                    mn[xi]  = f(X[xi], xt, l, v1.get('mn'), v2.get('mn'))
                    for ei in range(len(E)) : 
                        VA[ei][xi] = f(X[xi], xt, l, va1[ei], va2[ei])
    return T, B, ni, nn, nt, Xi, mi, mn, VA


