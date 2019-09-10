#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 11:14:36 2019

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt

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

# We define the constants of our problem
nt    = 0.35  # [atom/cm^-3] Density, here WNM
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
    
# We calculate the SNR Shock velocity
u_sh = np.empty(len(r_new))
for ii in range(1, len(r_new)) : 
    u_sh[ii] = (r_new[ii] - r_new[ii-1])/(t_new[ii] - t_new[ii-1])
u_sh[0] = u_sh[1] - (u_sh[2] - u_sh[1])
###############################################################################



# We calculate Emax(t) 

def f1(x, const) :
    a = const[0]
    c = const[1]
    return x*np.log(x/a) - c

def df1dx(x, const) : 
    a = const[0]
    return np.log(x/a)


def f2(x, const) : 
    a = const[0]
    b = const[1]
    c = const[2]
#    print ((x/a)**2)
    return x*((x/a)**b - 1.) - c

def df2dx(x, const) : 
    a = const[0]
    b = const[1]
    return b*(x/a)**b + (x/a)**b - 1.


def NewtonRaphson(f, df, x0, eps, const) : 
    exp = np.inf
    niter = 0
    x = x0 
    while (exp > eps) : 
        xold = x 
        x = x - f(x, const)/df(x, const)
        exp = abs(xold-x)/x
        #print (x/cst.GeV)
        niter += 1

    return x, niter


gamma = 2.2
beta = gamma - 2
Emin = 0.1*cst.GeV

Emax = np.empty(len(t_new))
eps = 1e-4
x0  = 10.*cst.GeV

if (beta != 0.) : 
    for ii in range(len(t_new)) : 
        a = Emin
        b = beta
        c = (beta/(1+beta))*cst.e*np.sqrt(4*np.pi*nt*cst.mp)/(10.*cst.c)*xhi_cr*u_sh[ii]**2*r_new[ii] 
        Emax[ii], niter = NewtonRaphson(f2, df2dx, x0, eps, [a, b, c])
#        print (niter)
if (beta == 0.) : 
    for ii in range(len(t_new)) : 
        a = Emin
        b = cst.e*np.sqrt(4*np.pi*nt*cst.mp)/(10.*cst.c)*xhi_cr*u_sh[ii]**2*r_new[ii]
        Emax[ii], niter = NewtonRaphson(f1, df1dx, x0, eps, [a, b, c])


Emin = 0.001*cst.GeV
EMAX = max(Emax)
delta = 3.

#tSed = 1e3*cst.yr*E51**(-0.5)*(Mej/10.)**(5./6)*(nt)**(-1./3) # [Celli et al. 2019]
tSed = tfree # [Truelove & McKee 1997]

def Gettesc(E, delta) : 
    return tSed*((E**2/cst.c**2 - cst.mp**2*cst.c**2)/(EMAX**2/cst.c**2 - cst.mp**2*cst.c**2))**(-1./(2.*delta))
    
Ecr = np.logspace(np.log10(0.1*cst.GeV), np.log10(100.*cst.TeV), 100)
tesc = np.empty(len(Ecr))
for ii in range(len(Ecr)) : 
    tesc[ii] = Gettesc(Ecr[ii], delta)




###############################################################################
# Model figure                                                                #
###############################################################################
plt.figure(figsize=(10,6))
for ii in range(len(t)) : 
    plt.loglog(t[ii]/cst.kyr, R[ii]/cst.pc, c="black", marker='o')
plt.loglog(t_new/cst.kyr, r_new/cst.pc, c="blue", lw=2, label="Spline model")

plt.legend()
plt.xlabel("Time [kyr]")
plt.ylabel("$R$ [pc]")


plt.figure(figsize=(10,6))
plt.loglog(t_new/cst.kyr, u_sh/(100*1e3), c="blue", lw=2, label="Shock Velocity")
plt.legend()
plt.ylabel("$V_\\mathrm{sh}$ [km/s]")
plt.xlabel("Time [kyr]")


plt.figure(figsize=(10,6))
plt.loglog(t_new/cst.kyr, Emax/cst.GeV, c="blue", lw=2, label="$\\Gamma = "+str(round(gamma,1))+" $")
plt.axvline(tSed/cst.kyr, c="black")
plt.legend()
plt.ylabel("$E_\\mathrm{max,0}(t)$ [GeV]")
plt.xlabel("Time [kyr]")

plt.figure(figsize=(10,6))
plt.loglog(Ecr/cst.GeV, tesc/cst.kyr, c="blue", lw=2, label="Escape Time")
plt.legend()
plt.ylabel("$t_\\mathrm{esc}$ [kyr]")
plt.xlabel("E [GeV]")

