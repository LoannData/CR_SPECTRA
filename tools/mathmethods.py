#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:04:53 2018

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import math
import scipy.special as sp 

import constants as cst 


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


# Routine d'integration de simpson avec un pas en log.
def g(f,t) :
    if (math.isnan(t) == True) : 
        return 0.0
    else: 
        return np.exp(t)*f(np.exp(t))
    
def simpson_log(f,a,b,N) : 
    a = np.log(a)
    b = np.log(b)
    S = 0. 
#    N = 1e3 # Nombre d'intervales de largeur h
    h = (b-a)/N 
    t = a
# Première partie    
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*g(f,t1)
    S2 = h*(17/48.)*g(f,t2)
    S3 = h*(59/48.)*g(f,t3)
    S4 = h*(43/48.)*g(f,t4)
    S5 = h*(49/48.)*g(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
# Boucle pour les termes suivants    
    while (t <= b-4*h) : 
        Si = h*g(f,t)
        t += h
        S += Si   
# Dernière partie pour les derniers termes 
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*g(f,t1)
    S2 = h*(49/48.)*g(f,t2)
    S3 = h*(43/48.)*g(f,t3)
    S4 = h*(59/48.)*g(f,t4)
    S5 = h*(17/48.)*g(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
        
    return S

# Routine d'integration de simpson avec un pas linéaire
def glin(f,t) : 
    return f(t)
    
def simpson_lin(f,a,b,N) : 
#    a = np.log(a)
#    b = np.log(b)
    S = 0. 
#    N = 1e3 # Nombre d'intervales de largeur h
    h = (b-a)/N 
    t = a
# Première partie    
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*glin(f,t1)
    S2 = h*(17/48.)*glin(f,t2)
    S3 = h*(59/48.)*glin(f,t3)
    S4 = h*(43/48.)*glin(f,t4)
    S5 = h*(49/48.)*glin(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
# Boucle pour les termes suivants    
    while (t <= b-4*h) : 
        Si = h*glin(f,t)
        t += h
        S += Si   
# Dernière partie pour les derniers termes 
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*glin(f,t1)
    S2 = h*(49/48.)*glin(f,t2)
    S3 = h*(43/48.)*glin(f,t3)
    S4 = h*(59/48.)*glin(f,t4)
    S5 = h*(17/48.)*glin(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
        
    return S


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

def shape(X, Amp, Xmin = 0., Xmax = 1., sig_Xmin = 0.1, sig_Xmax = 0.2) : 
    if (X - Xmin == 0. or Xmax - X == 0.) : 
        return Amp
    return 0.5*Amp*(sp.erf((X - Xmin)/sig_Xmin) + sp.erf((Xmax - X)/sig_Xmax))

def multishape(X, Amp, Xmin, Xmax, sig) : 
    sig_Xmin = [0] + sig
    sig_Xmax = sig + [0]
    S = 0. 
    for i in range(len(Amp)) : 
        S += shape(X, Amp[i], Xmin[i], Xmax[i], sig_Xmin[i], sig_Xmax[i])
    return S

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
    gamma_in = np.zeros((len(E), len(X)))
    gamma_lz = np.zeros((len(E), len(X)))
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
                        VA[ei][xi]         = phases[pi][2][ei]
                        gamma_in[ei][xi]   = phases[pi][3][ei]
                        gamma_lz[ei][xi]   = phases[pi][4][ei]
    if (len(phases) > 1) : 
        print ("Multiphase setup : Calculation of the smoothing parameters")
        for ei in range(len(E)) : 
            print ("Calculation of the smoothing parameters - Energy index ",ei," over ",len(E))
            Xmin   = []
            Xmax   = []
            Amp_Va       = []
            Amp_gamma_in = []
            Amp_gamma_lz = []
            for pi in range(len(phases)) :
                Xmin.append(phases[pi][1].get("Xmin"))
                Xmax.append(phases[pi][1].get("Xmax"))
                Amp_Va.append(phases[pi][2][ei])
                Amp_gamma_in.append(phases[pi][3][ei])
                Amp_gamma_lz.append(phases[pi][4][ei])
            for xi in range(len(X)) : 
                VA[ei][xi] = multishape(X[xi], Amp_Va, Xmin, Xmax, smooth_width)
                gamma_in[ei][xi] = multishape(X[xi], Amp_gamma_in, Xmin, Xmax, smooth_width)
                gamma_lz[ei][xi] = multishape(X[xi], Amp_gamma_lz, Xmin, Xmax, smooth_width)
            
        
        Xmin   = []
        Xmax   = []
        Amp_T  = []
        Amp_B  = []
        Amp_ni = []
        Amp_nn = []
        Amp_nt = []
        Amp_Xi = []
        Amp_mi = []
        Amp_mn = []
        for pi in range(len(phases)) :
            Xmin.append(phases[pi][1].get("Xmin"))
            Xmax.append(phases[pi][1].get("Xmax"))
            Amp_T.append(phases[pi][0].get("T"))
            Amp_B.append(phases[pi][0].get("B"))
            Amp_ni.append(phases[pi][0].get("ni"))
            Amp_nn.append(phases[pi][0].get("nn"))
            Amp_nt.append(phases[pi][0].get("nt"))
            Amp_Xi.append(phases[pi][0].get("X"))
            Amp_mi.append(phases[pi][0].get("mi"))
            Amp_mn.append(phases[pi][0].get("mn"))
        for xi in range(len(X)) : 
            T[xi] = multishape(X[xi], Amp_T, Xmin, Xmax, smooth_width)
            B[xi] = multishape(X[xi], Amp_B, Xmin, Xmax, smooth_width)
            ni[xi] = multishape(X[xi], Amp_ni, Xmin, Xmax, smooth_width)
            nn[xi] = multishape(X[xi], Amp_nn, Xmin, Xmax, smooth_width)
            nt[xi] = multishape(X[xi], Amp_nt, Xmin, Xmax, smooth_width)
            Xi[xi] = multishape(X[xi], Amp_Xi, Xmin, Xmax, smooth_width)
            mi[xi] = multishape(X[xi], Amp_mi, Xmin, Xmax, smooth_width)
            mn[xi] = multishape(X[xi], Amp_mn, Xmin, Xmax, smooth_width)

    return T, B, ni, nn, nt, Xi, mi, mn, VA, gamma_in, gamma_lz

"""
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
    gamma_in = np.zeros((len(E), len(X)))
    gamma_lz = np.zeros((len(E), len(X)))
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
                        VA[ei][xi]         = phases[pi][2][ei]
                        gamma_in[ei][xi]   = phases[pi][3][ei]
                        gamma_lz[ei][xi]   = phases[pi][4][ei]
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
                        gamma_in1 = phases[pos][3]
                        gamma_in2 = phases[pos][3]
                        gamma_lz1 = phases[pos][4]
                        gamma_lz2 = phases[pos][4]
                    else : 
                        v1 = phases[pos-1][0]
                        v2 = phases[pos][0]
                        va1 = phases[pos-1][2]
                        va2 = phases[pos][2]
                        gamma_in1 = phases[pos-1][3]
                        gamma_in2 = phases[pos][3]
                        gamma_lz1 = phases[pos-1][4]
                        gamma_lz2 = phases[pos][4]
                    T[xi]  = f(X[xi], xt, l, v1.get('T'), v2.get('T'))
                    B[xi]  = f(X[xi], xt, l, v1.get('B'), v2.get('B'))
                    ni[xi]  = f(X[xi], xt, l, v1.get('ni'), v2.get('ni'))
                    nn[xi]  = f(X[xi], xt, l, v1.get('nn'), v2.get('nn'))
                    nt[xi] = f(X[xi], xt, l, v1.get('nt'), v2.get('nt'))
                    Xi[xi]  = f(X[xi], xt, l, v1.get('X'), v2.get('X'))
                    mi[xi]  = f(X[xi], xt, l, v1.get('mi'), v2.get('mi'))
                    mn[xi]  = f(X[xi], xt, l, v1.get('mn'), v2.get('mn'))
                    for ei in range(len(E)) : 
                        VA[ei][xi]       = f(X[xi], xt, l, va1[ei], va2[ei])
                        gamma_in[ei][xi] = f(X[xi], xt, l, gamma_in1[ei], gamma_in2[ei])
                        gamma_lz[ei][xi] = f(X[xi], xt, l, gamma_lz1[ei], gamma_lz2[ei])
                if (X[xi] >= phases[pos][1].get('Xmin') + 0.5*(phases[pos][1].get('Xmax') - phases[pos][1].get('Xmin')) and X[xi] <= phases[pos][1].get('Xmax')) :
                    xt = phases[pos][1].get('Xmax')
                    l  = smooth_width
                    if (pos == len(phases)-1) :
                        v1 = phases[pos][0]
                        v2 = phases[pos][0]
                        va1 = phases[pos][2]
                        va2 = phases[pos][2]
                        gamma_in1 = phases[pos][3]
                        gamma_in2 = phases[pos][3]
                        gamma_lz1 = phases[pos][4]
                        gamma_lz2 = phases[pos][4]
                    else : 
                        v1 = phases[pos][0]
                        v2 = phases[pos+1][0]
                        va1 = phases[pos][2]
                        va2 = phases[pos+1][2]
                        gamma_in1 = phases[pos][3]
                        gamma_in2 = phases[pos+1][3]
                        gamma_lz1 = phases[pos][4]
                        gamma_lz2 = phases[pos+1][4]
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
                        gamma_in[ei][xi] = f(X[xi], xt, l, gamma_in1[ei], gamma_in2[ei])
                        gamma_lz[ei][xi] = f(X[xi], xt, l, gamma_lz1[ei], gamma_lz2[ei])
    return T, B, ni, nn, nt, Xi, mi, mn, VA, gamma_in, gamma_lz
"""

