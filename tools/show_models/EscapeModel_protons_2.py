#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:40:12 2020

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys 
sys.path.append("../")
import constants as cst
import phases_collection as ism 

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
###############################################################################


import sys 
sys.path.append('../tools/')


def getSNR(phase, size = 100) : 

    # We define the constants of our problem
    nt    = phase.get("ni") + phase.get("nn")
    xhi_m = 1 # Metallicity (1 for solar abundances)
    
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
    logt_new = np.linspace(logt[0], logt[-1], size)
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
    
    
    
    return {"R_SNR"   :r_new, 
            "t_SNR"   :t_new, 
            "u_sh"    :u_sh, 
            "t_SED"   :tfree, 
            "t_PDS"   :tPDS,
            "t_MCS"   :tMCS,
            "t_merge" :tmerge,
            "t_max"   : tmax, 
            "r_SED"   :R_free,
            "r_PDS"   : R_PDS,
            "r_MCS"   : R_MCS,
            "r_merge" : R_merge}


def Gettesc(E, delta, tSed, EMAX, Emin = 0.1*cst.GeV) : 
    # E = E + 1.*cst.GeV
    # E = E - 1*cst.GeV
    if (E < EMAX and E > Emin) : 
        # return tSed*((E**2/cst.c**2 - cst.mp**2*cst.c**2)/(EMAX**2/cst.c**2 - cst.mp**2*cst.c**2))**(-1./(2.*delta))
        return tSed*((E**2/cst.c**2)/(EMAX**2/cst.c**2))**(-1./(2.*delta))
    else : 
        return np.NaN


def getEmax(t_new, u_sh, r_new, phase, gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV) : 
    xhi_cr= 0.1
    nt    = phase.get("ni") + phase.get("nn")
    beta = gamma - 2
    

    Emax = np.empty(len(t_new))
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
    
    return Emax 


Ecr = np.logspace(np.log10(0.1*cst.GeV), np.log10(100.*cst.TeV), 1000)


SNR_HII  = getSNR(ism.HII, size = 100)
emax_HII = getEmax(SNR_HII.get("t_SNR"), SNR_HII.get("u_sh"), SNR_HII.get("R_SNR"), ism.HII, 
                   gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV)

SNR_WIM  = getSNR(ism.WIM, size = 100)
emax_WIM = getEmax(SNR_WIM.get("t_SNR"), SNR_WIM.get("u_sh"), SNR_WIM.get("R_SNR"), ism.WIM, 
                   gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV)


SNR_WNM  = getSNR(ism.WNM, size = 100)
emax_WNM = getEmax(SNR_WNM.get("t_SNR"), SNR_WNM.get("u_sh"), SNR_WNM.get("R_SNR"), ism.WNM, 
                   gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV)

SNR_CNM  = getSNR(ism.CNM, size = 100)
emax_CNM = getEmax(SNR_CNM.get("t_SNR"), SNR_CNM.get("u_sh"), SNR_CNM.get("R_SNR"), ism.CNM, 
                   gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV)

SNR_DiM  = getSNR(ism.DiM, size = 100)
emax_DiM = getEmax(SNR_DiM.get("t_SNR"), SNR_DiM.get("u_sh"), SNR_DiM.get("R_SNR"), ism.DiM, 
                   gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV)


delta = 2. 
tesc_HII = np.empty(len(Ecr))
tesc_WIM = np.empty(len(Ecr))
tesc_WNM = np.empty(len(Ecr))
tesc_CNM = np.empty(len(Ecr))
tesc_DiM = np.empty(len(Ecr))

tesc_HII_3 = np.empty(len(Ecr))
tesc_WIM_3 = np.empty(len(Ecr))
tesc_WNM_3 = np.empty(len(Ecr))
tesc_CNM_3 = np.empty(len(Ecr))
tesc_DiM_3 = np.empty(len(Ecr))

for ii in range(len(Ecr)) : 
    tesc_HII[ii] = Gettesc(Ecr[ii], delta, SNR_HII.get("t_SED"), max(emax_HII))
    tesc_WIM[ii] = Gettesc(Ecr[ii], delta, SNR_WIM.get("t_SED"), max(emax_WIM))
    tesc_WNM[ii] = Gettesc(Ecr[ii], delta, SNR_WNM.get("t_SED"), max(emax_WNM))
    tesc_CNM[ii] = Gettesc(Ecr[ii], delta, SNR_CNM.get("t_SED"), max(emax_CNM))
    tesc_DiM[ii] = Gettesc(Ecr[ii], delta, SNR_DiM.get("t_SED"), max(emax_DiM))
    
    tesc_HII_3[ii] = Gettesc(Ecr[ii], 3., SNR_HII.get("t_SED"), max(emax_HII))
    tesc_WIM_3[ii] = Gettesc(Ecr[ii], 3., SNR_WIM.get("t_SED"), max(emax_WIM))
    tesc_WNM_3[ii] = Gettesc(Ecr[ii], 3., SNR_WNM.get("t_SED"), max(emax_WNM))
    tesc_CNM_3[ii] = Gettesc(Ecr[ii], 3., SNR_CNM.get("t_SED"), max(emax_CNM))
    tesc_DiM_3[ii] = Gettesc(Ecr[ii], 3., SNR_DiM.get("t_SED"), max(emax_DiM))


size_x = 4
size_y = 3.5
sub_x  = 2
sub_y  = 1
fig = plt.figure(figsize=(size_x*sub_x,size_y*sub_y))

gs = gridspec.GridSpec(ncols= sub_x, nrows = sub_y, figure = fig )
gs.update(wspace=0.05, hspace=0.05) # set the spacing between axes.


ax0 = fig.add_subplot(gs[0])

ax0.loglog(SNR_HII.get("t_SNR")/cst.kyr, emax_HII/cst.GeV, c="red", label="HII")
ax0.loglog(SNR_WIM.get("t_SNR")/cst.kyr, emax_WIM/cst.GeV, c="orange", label ="WIM")
ax0.loglog(SNR_WNM.get("t_SNR")/cst.kyr, emax_WNM/cst.GeV, c="green", label="WNM")
ax0.loglog(SNR_CNM.get("t_SNR")/cst.kyr, emax_CNM/cst.GeV, c="deepskyblue", label="CNM")
ax0.loglog(SNR_DiM.get("t_SNR")/cst.kyr, emax_DiM/cst.GeV, c="blue", label="DiM")

ax0.legend(loc="upper right", ncol = 5, bbox_to_anchor=(1.8, 1.15))

ax0.set_ylabel("$E_{\\mathrm{max},0}$ [GeV]")
ax0.set_xlabel("$t$ [kyr]")


ax1 = fig.add_subplot(gs[1])

ax1.loglog(Ecr/cst.GeV, tesc_HII/cst.kyr, c="red", ls ="-")
ax1.loglog(Ecr/cst.GeV, tesc_WIM/cst.kyr, c="orange", ls ="-")
ax1.loglog(Ecr/cst.GeV, tesc_WNM/cst.kyr, c="green", ls ="-")
ax1.loglog(Ecr/cst.GeV, tesc_CNM/cst.kyr, c="deepskyblue", ls ="-")
ax1.loglog(Ecr/cst.GeV, tesc_DiM/cst.kyr, c="blue", ls ="-")

ax1.loglog(Ecr/cst.GeV, tesc_HII_3/cst.kyr, c="red", ls ="--")
ax1.loglog(Ecr/cst.GeV, tesc_WIM_3/cst.kyr, c="orange", ls ="--")
ax1.loglog(Ecr/cst.GeV, tesc_WNM_3/cst.kyr, c="green", ls ="--")
ax1.loglog(Ecr/cst.GeV, tesc_CNM_3/cst.kyr, c="deepskyblue", ls ="--")
ax1.loglog(Ecr/cst.GeV, tesc_DiM_3/cst.kyr, c="blue", ls ="--")

ax1.plot([],[],ls='-',c="black",label="$\\delta = 2$")
ax1.plot([],[],ls='--',c="black",label="$\\delta = 3$")

ax1.legend(loc="upper right")

ax1.yaxis.tick_right()
ax1.yaxis.set_label_position("right")
ax1.set_ylabel("Escape time $t_\\mathrm{esc}$ [kyr]")
ax1.set_xlabel("$E$ [GeV]")
ax1.set_xlim(1e-1, 1e5)
ax1.set_ylim(4e-2, 1e2)


fig.tight_layout()
fig.savefig("escape_model.pdf",pad=-10)










