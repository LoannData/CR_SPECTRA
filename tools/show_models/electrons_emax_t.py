#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:36:52 2020

@author: lbrahimi
"""
import numpy as np 
import matplotlib.pyplot as plt

import sys 
sys.path.append("../")
import constants as cst

###############################################################################
# FUNCTIONS IN ORDER TO CREATE A GOOD SPLINE !!!                              #
###############################################################################
# TriDiagonal matrix inversion function
def InverseTrigonalMatrix(T) : 
    n = len(T)
    
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    for i in range(n) : 
        b[i] = T[i][i]
        if (i < n-1) : c[i] = T[i][i+1]
        if (i > 0) : a[i] = T[i][i-1]
    
    theta_mm  = 0
    theta_m   = 1
    theta = np.zeros(n)
    theta[0] = b[0]*theta_m #- a[0]*c[0]*theta_mm
    theta[1] = b[1]*theta[0] - a[1]*c[0]*theta_m
    for i in range(2, n, 1) : 
        theta[i] = b[i]*theta[i-1] - a[i]*c[i-1]*theta[i-2]
    
    phi_pp = 0
    phi_p  = 1
    phi = np.zeros(n+2)
    phi[n] = phi_p
    phi[n+1] = phi_pp
    phi[n-1] = b[n-1]*phi_p #- c[n-1]*a[n]*phi_pp
    phi[n-2] = b[n-2]*phi[n-1] - c[n-2]*a[n-1]*phi_p
    for i in range(n-3, -1, -1) : 
        phi[i] = b[i]*phi[i+1] - c[i]*a[i+1]*phi[i+2]
        
    
    Tinv = np.empty((n, n))
    for i in range(n) : 
        for j in range(n) : 
    #        print (i,j)
            if (i < j) : 
                p = 1.
                for ei in range(i, j, 1) : 
                    p = p*c[ei]
                if (i-1 == -1) : 
                    loc_theta = theta_m
                if (i-1 > -1) : 
                    loc_theta = theta[i-1]
                
                Tij = (-1)**(i+j)*p*loc_theta*phi[j+1]/theta[n-1]
            if (i == j) : 
                if (i-1 == -1) : 
                    loc_theta = theta_m
                if (i-1 > -1) : 
                    loc_theta = theta[i-1]
                Tij = loc_theta*phi[j+1]/theta[n-1]
            if (i > j) : 
                p = 1.
                for ei in range(j+1, i+1, 1) : 
                    p = p*a[ei]
                if (j-1 == -1) : 
                    loc_theta = theta_m
                if (j-1 > -1) : 
                    loc_theta = theta[j-1]
                Tij = (-1)**(i+j)*p*loc_theta*phi[i+1]/theta[n-1]
            Tinv[i, j] = Tij
    
    return Tinv 


def ProductMatrix(A, B) : 
    A_l = len(A)
    A_c = len(A[0])
    
    B_l = len(B)
    B_c = len(B[0])
    
    C_l = A_l
    C_c = B_c
    
    C = np.zeros((C_l, C_c))
    for i in range(C_l) : 
        for j in range(C_c) : 
            s = 0.
            for k in range(A_c) : 
                s += A[i][k]*B[k][j]
            C[i][j]=  s 

    return C


def InterpolatingSpline(X, Y) : 

    N = len(X)
#    print (N)

    h = []
    for i in range(N-1) : 
        h.append([X[i+1] - X[i]])
    
    lowF  = 0.
    highF = 0.
    F = [[lowF]]
    for i in range(1, N-1) : 
        F.append([(Y[i+1] - Y[i])/h[i] - (Y[i] - Y[i-1])/h[i-1]])
    F.append([highF])

    R = np.zeros((N, N))
    R[0][0]     = 1.
    R[N-1][N-1] = 1.
    
    for i in range(1, N-1) : 
        R[i][i]   = (h[i-1][0] + h[i][0])/3.
        R[i][i+1] = h[i][0]/6.
        R[i][i-1] = h[i-1][0]/6.
    

    R_inv = InverseTrigonalMatrix(R)

    M = ProductMatrix(R_inv, F)
    
    Ci  = np.empty(N-1)
    Cip = np.empty(N-1)
    
    for i in range(N-1) : 
        Ci[i]  = (Y[i+1] - Y[i])/h[i][0] - h[i][0]*(M[i+1][0] - M[i][0])/6.
        Cip[i] = Y[i] - M[i][0]*h[i][0]**2/6.
    
    
    def f(x) : 
        for k in range(0, N-1) : 
            if (x >= X[k] and x <= X[k+1]) : 
                loc_a = M[k][0]*(X[k+1] - x)**3/(6.*h[k][0])
                loc_b = M[k+1][0]*(x - X[k])**3/(6.*h[k][0])
                loc_c = Ci[k]*(x-X[k]) + Cip[k]
        return loc_a + loc_b + loc_c 
    
    return f
###############################################################################
    
import sys 
sys.path.append('../tools/')
import constants as cst


def Rsh(nt) : 

    xhi_m = 1 # Metallicity (1 for solar abundances)
    xhi_cr= 0.1
    E51   = 1 # [10^51 erg] Total energy released by the SNR
    Mej   = 1 # [Msum] Total mass released by the SNR
    C06   = 1 # 
    beta  = 1 #
    phi_c = 1 # Ratio of the thermal conductivity to the Spitzer (1962) value
    vej8  = 10.*(E51/Mej)**(0.5)
    ###############################################################################
    # FUNCTIONS IN ORDER TO MAKE OUR SNR EXPAND IN THE ISM                        #
    ###############################################################################
    # We define the characteristic times of our problem
    tini   = 1e-4*cst.kyr # [s]
    tfree  = 0.3*E51**(-0.5)*Mej*nt**(-1./3)*cst.kyr # [s]
    tPDS   = np.exp(-1.)*3.61e4*E51**(3./14)/(xhi_m**(5./14)*nt**(4./7))*cst.yr # [s]
    tMCS   = min(61*vej8**3/(xhi_m**(9./14)*nt**(3./7)*E51**(3./14)), 476./(xhi_m*phi_c)**(9./14))*tPDS # [s]
    tmerge = 153.*(E51**(1./14)*nt**(1./7)*xhi_m**(3./14)/(beta*C06))**(10./7)*tPDS # [s]
    tmax = min(tMCS, tmerge) # [s]
    
    # We define the characteristic radii of our problem
    R_free    = 5.0*(E51/nt)**(1./5)*(1 - (0.05*Mej**(5./6))/(E51**0.5*nt**(1./3)*(tfree/cst.kyr)))**(2./5)*(tfree/cst.kyr)**(2./5)*cst.pc # [cm]
    R_ini     = R_free*(tini/tfree)**(1.)
    R_PDS     = 5.0*(E51/nt)**(1./5)*(1 - (0.05*Mej**(5./6))/(E51**0.5*nt**(1./3)*(tPDS/cst.kyr)))**(2./5)*(tPDS/cst.kyr)**(2./5)*cst.pc # [cm] 
    R_MCS     = R_PDS*(tMCS/tPDS)**(3./10)
    R_merge   = R_MCS*(tmerge/tMCS)**(1./4)
    if (tMCS < tmerge) : 
        t = np.array([tini, tfree, tPDS, tMCS, tmerge])
        R = np.array([R_ini, R_free, R_PDS, R_MCS, R_merge])
    if (tMCS >= tmerge) : 
        t = np.array([tini, tfree, tPDS, tmerge])
        R = np.array([R_ini, R_free, R_PDS, R_merge])
    
    # We pass in the loglog space 
    logt = np.empty(len(t))
    logR = np.empty(len(R))
    for ii in range(len(t)) : 
        logt[ii] = np.log10(t[ii])
        logR[ii] = np.log10(R[ii])
    # We apply the spline interpolation
    f_SNR = InterpolatingSpline(logt, logR)
    logt_new = np.linspace(logt[0], logt[-1], 100)
    logr_new = np.empty(len(logt_new))
    for ii in range(len(logt_new)) : 
        logr_new[ii] = f_SNR(logt_new[ii])
    # We come back in the linear space 
    t_new = np.empty(len(logt_new))
    r_new = np.empty(len(logr_new))
    for ii in range(len(t_new)) : 
        t_new[ii] = 10**(logt_new[ii])
        r_new[ii] = 10**(logr_new[ii])
        
    return t, R, t_new, r_new



###############################################################################
# ELECTRON ESCAPE MODEL (from Ohira et al. (2012))                            # 
###############################################################################


# Some important functions 

def B(t, ts, tB, alpha_B, BISM, Bfree) :
    if (t <= ts) : 
        return Bfree
    if (t > ts and t <= tB) : 
        return Bfree*(t/ts)**(-alpha_B)
    if (t > tB) : 
        return BISM


def eta_g(t, ts, tB, alpha, alpha_B, eta_free) : 
    if (t <= ts) : 
        return eta_free
    if (t > ts and t <= tB) : 
        return eta_free*(t/ts)**(alpha - alpha_B - 1./5)
    if (t > tB) : 
        return eta_free*(tB/ts)**(-alpha_B)*(t/ts)**(alpha - 1./5)

def Em_age(t, ts, tB, alpha_B, alpha, BISM, Bfree, eta_acc, eta_free, Rs) : 
    C = 3*cst.e*B(t, ts, tB, alpha_B, BISM, Bfree)*Rs**2/(eta_acc*eta_g(t, ts, tB, alpha, alpha_B, eta_free)*cst.c*ts)
    if (t <= ts) : 
        return C*(t/ts)
    if (t > ts) : 
        return C*(t/ts)**(-1./5)

def Em_cool(t, ts, tB, alpha_B, alpha, BISM, Bfree, Ems) : 
    if (t <= ts) : 
        return Ems
    if (t > ts and t <= tB) : 
        return Ems*(t/ts)**((2*alpha_B - alpha - 1)/2.)
    if (t > tB) : 
        return (tB/ts)**(alpha_B)*(t/ts)**(-(alpha + 1)/2.)

def Em_esc(t, ts, tB, alpha_B, alpha, BISM, Bfree, eta_acc, eta_free, eta_esc, Rs) : 
    return np.sqrt(eta_esc*eta_acc)*Em_age(t, ts, tB, alpha_B, alpha, BISM, Bfree, eta_acc, eta_free, Rs)




def Emax_electrons(time) : 

# time = 12.5*cst.kyr 

    E51   = 1 # [10^51 erg] Total energy released by the SNR
    Mej   = 1 # [Msum] Total mass released by the SNR
    nt    = 0.35 # [cm^{-3}] Total ISM density 
    B0    = 6e-6 # [G] Mean magnetic field density in the ISM
    
    eta_acc   = 10  # Numerical factor depending on the shock compression ratio 
    eta_free = 1   # Gyrofactor during the free expansion phase
    eta_esc   = 0.1 # Numerical factor for the escape time 
    Eknee     = 10**(15.5)*cst.eV
    
    # alpha_B coefficient 
    alpha = 2.6
    # alpha_B   =  alpha - 1./5
    alpha_B   =  9./10
    # alpha_B   =  3./5
    
    
    # We define the time corresponding to the beginning of the Sedov phase 
    ts  = 0.3*E51**(-0.5)*Mej*nt**(-1./3)*cst.kyr # [s]
    Rs = 5.0*(E51/nt)**(1./5)*(1 - (0.05*Mej**(5./6))/(E51**0.5*nt**(1./3)*(ts/cst.kyr)))**(2./5)*(ts/cst.kyr)**(2./5)*cst.pc
    
    # Amplified magnetic field during the free expansion phase 
    Bfree = eta_free*eta_acc*cst.c*ts*Eknee/(3*cst.e*Rs**2) 
    
    tB = ts*(Bfree/B0)**(1./alpha_B)
    
    Em_s = 9*cst.me**2*cst.c**(5./2)*Rs**2/(8*eta_free*eta_acc*cst.e*ts**(3./2)*Eknee**(1./2))
    
    
    EM_age  = Em_age(time, ts, tB, alpha_B, alpha, B0, Bfree, eta_acc, eta_free, Rs)
    EM_cool = Em_cool(time, ts, tB, alpha_B, alpha, B0, Bfree, Em_s)
    EM_esc  = Em_esc(time, ts, tB, alpha_B, alpha, B0, Bfree, eta_acc, eta_free, eta_esc, Rs)
    
    Em_e = min(min(EM_age, EM_cool), min(EM_cool, EM_esc))

    return Em_e 




def escape_time(E, tSed, EM, delta) : 
    return tSed*((E**2/cst.c**2 - cst.me**2*cst.c**2)/(EM**2/cst.c**2 - cst.me**2*cst.c**2))**(-1/(2*delta))



# time = np.logspace(np.log10(1e-3*cst.kyr), np.log10(1e3*cst.kyr), num=100)
# Em_e = np.zeros(len(time)) 
# for ii in range(len(time)) : 
#     Em_e[ii] = Emax_electrons(time[ii])
# plt.figure()
# plt.loglog(time/cst.kyr, Em_e/cst.GeV)
    
# alpha_B coefficient 
alpha = 2.6
# alpha_B   =  alpha - 1./5
alpha_B   =  9./10
# alpha_B   =  3./5



xhi_m = 1 # Metallicity (1 for solar abundances)
xhi_cr= 0.1
E51   = 1 # [10^51 erg] Total energy released by the SNR
Mej   = 1 # [Msum] Total mass released by the SNR
C06   = 1 # 
beta  = 1 #
phi_c = 1 # Ratio of the thermal conductivity to the Spitzer (1962) value
vej8  = 10.*(E51/Mej)**(0.5)

nt = 0.35
###############################################################################
# FUNCTIONS IN ORDER TO MAKE OUR SNR EXPAND IN THE ISM                        #
###############################################################################
# We define the characteristic times of our problem
tsed  = 0.3*E51**(-0.5)*Mej*nt**(-1./3)*cst.kyr # [s]


EM = Emax_electrons(tsed)


E     = np.logspace(np.log10(1*cst.GeV), np.log10(100*cst.TeV), num=100)
tesc  = np.zeros(len(E))

for ii in range(len(E)) :
    tesc[ii] = escape_time(E[ii], tsed, EM, 0.25*(alpha + 2*alpha_B - 1))


plt.figure()
plt.loglog(E/cst.GeV, tesc/cst.kyr)




































