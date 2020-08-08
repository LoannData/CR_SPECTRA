#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:12:16 2020

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
import damping as dp 

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
            "r_SED"   : R_free,
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


def IonNeutral_Damping(k, medium_props, nu_n = 0, theta = 0) : 
    T  = medium_props.get('T')
    B  = medium_props.get('B')
    mn = medium_props.get('mn')
    mi = medium_props.get('mi') 
    X  = medium_props.get('X')
    ni = medium_props.get('ni')
    nn = medium_props.get('nn')    
    
    
    # Some relations 
    # rl = (E)/(cst.e*B)
    # k = 1/rl
    chi = (mn/mi)*(X**(-1.)-1.)
    
    VAi = B/np.sqrt(4*np.pi*mi*ni)
    if (T <= 50) : 
        nu_in = 2*nn*8.4e-9*(50/1e4)**0.4
    if (T > 50) : 
        nu_in = 2*nn*8.4e-9*(T/1e4)**0.4
    nu_ni = chi**(-1.)*nu_in
    
    kz = k*np.cos(theta)
    
    a = 1 
    b = (1 + chi)*nu_ni
    c = kz**2*VAi**2
    d = nu_ni*kz**2*VAi**2
    
    roots = mt.Cubic3(a, b, c, d)
    
    wR = abs(roots[2].imag)
    wI = abs(roots[2].real)
    cA = wR/k 
    
    return dict(wr=wR, wi=-wI, VA=cA)


def Duu_Alfven_Slab_Linear_Undamped(mu, E, medium_props, mass = cst.mp, kmin = 1e-20, q = 1.5, I = 1e-4) : 
    """
    See. Schlickeiser (2002, Chap 13.1.3.1, p.318)

    Parameters
    ----------
    mu : TYPE
        DESCRIPTION.
    medium_props : TYPE
        DESCRIPTION.
    particles_props : TYPE
        DESCRIPTION.
    kmin : TYPE, optional
        DESCRIPTION. The default is 1e-20.
    q : TYPE, optional
        DESCRIPTION. The default is 1.5.

    Returns
    -------
    D : TYPE
        DESCRIPTION.

    """
    
    def g_s(k, kmin, q, medium_props, I = I) : 
        B0 = medium_props.get("B")
        W0 = B0**2/(8*np.pi)
        I  = I
        g_s0 = 2*(q - 1)*W0*I*kmin**(q - 1)
        # g_s0 = I*W0*kmin
        if (k >= 0) : 
            if (k >= kmin) : 
                return g_s0*k**(-q)
            else : 
                return 0. 
        if (k < 0) : 
            if (abs(k) >= kmin) : 
                return g_s0*abs(k)**(-q)
            else : 
                return 0.
    
    
    m = mass 
    gamma = 1 + (E /(m*cst.c**2))
    v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
    p = gamma*m*v 
    Omega0 = cst.e*medium_props.get("B")/(m*cst.c)
    Omega  = Omega0/gamma 
    
    
    B0    = medium_props.get("B")
    mn    = medium_props.get("mn")
    nn    = medium_props.get("nn")
    mi    = medium_props.get("mi")
    ni    = medium_props.get("ni")
    
    # eps = B0/np.sqrt(4*np.pi*(mi*ni + mn*nn))/v
    eps = B0/np.sqrt(4*np.pi*(mi*ni))/v
    
    k_rp = (Omega/v)/(mu - eps)
    k_rm = (Omega/v)/(mu + eps)
    
    Cte1    = np.pi*Omega**2*(1 - mu**2)/(v*B0**2)
    Cte2_jp = (1 - eps*mu)**2/abs(mu - eps)
    Cte2_jm = (1 + eps*mu)**2/abs(mu + eps)
    Gp_jp      = g_s(k_rp, kmin, q, medium_props)  
    Gm_jp      = g_s(-k_rp, kmin, q, medium_props) 
    Gp_jm      = g_s(k_rm, kmin, q, medium_props)  
    Gm_jm      = g_s(-k_rm, kmin, q, medium_props) 
    
    Djp = Cte1*Cte2_jp*(Gp_jp + Gm_jp)
    Djm = Cte1*Cte2_jm*(Gp_jm + Gm_jm)
    
    D = Djp + Djm
    return D




def Kappa_zz(E, medium_props, mass = cst.mp, kmin = 1e-20, q = 5./3, I = 1e28) : 
    """
    

    Parameters
    ----------
    E : TYPE
        DESCRIPTION.
    medium_props : TYPE
        DESCRIPTION.
    mass : TYPE float, optional
        Mass of the diffusin particle. The defaut is m_proton
    kmin : TYPE float, optional
        Minimun length in cm^-1 for the turbulence spectra. The default is 1e-20.
    q : TYPE float, optional
        Spectral index of the Kolmogorov-like turbulence 
        spectrum. The default is 5./3.
    I : TYPE float, optional
        Diffusion coefficient normalization value for 1GeV particle and 5./3 spectrum. The default is 1e28.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    # I = (2.2345649450772952e28*1e6/I)
    
    
    
    m = mass
    gamma = 1 + (E /(m*cst.c**2))
    v = cst.c*np.sqrt(1 - (1/(E/(m*cst.c**2) + 1))**2)
    p = gamma*m*v 
    Omega0 = cst.e*medium_props.get("B")/(m*cst.c)
    Omega  = Omega0/gamma 
    
    lz_damping = dp.damping_lazarian_nopos(E, medium_props)
    # print (lz_damping)
    kmax = (lz_damping[1])**(-1)
    # kmax = 1e-21
    
    # print (kmax)
    
    I = I*kmax
    
    k_zz = 0.
    mu = np.linspace(-1, 1, 100)
    for jj in range(1, len(mu)-1) : 
        dmumu = 0.5*(mu[jj+1] - mu[jj-1])
        k_zz += dmumu*(1 - mu[jj]**2)**2/Duu_Alfven_Slab_Linear_Undamped(mu[jj], E, medium_props, mass = m, kmin = kmin, q = q, I = I)
    
    return v**2/8*k_zz



def kappa_zz_BC(E, d00, delta) : 
    return d00*(E/(10*cst.GeV))**(delta)



SNR_phase = ism.WNM
ISM_phase = [ism.WNM, ism.CNM, ism.DiM]
ISM_color = ["orange", "green", "blue"]
ISM_name  = ["WNM", "CNM", "DiM"]
delta = 4 



E = np.logspace(np.log10(0.1*cst.GeV), np.log10(100.*cst.TeV), 1000)

kzz_i = np.zeros(len(E)) 
kzz_d = np.zeros((len(ISM_phase),len(E)))
tesc_SNR = np.empty(len(E))
Resc     = np.empty(len(E))
SNR  = getSNR(SNR_phase, size = 1000)
emax_SNR = getEmax(SNR.get("t_SNR"), SNR.get("u_sh"), SNR.get("R_SNR"), SNR_phase, 
                   gamma = 2.2, Emin = 0.1*cst.GeV, eps = 1e-4, x0 = 10.*cst.GeV)

for ii in range(len(E)) : 
    tesc_SNR[ii] = Gettesc(E[ii], delta, SNR.get("t_SED"), max(emax_SNR))
    kzz_i[ii] =  kappa_zz_BC(E[ii], 1e28, 0.5)
    for ni in range(len(ISM_phase)) : 
        kzz_d[ni][ii] =  Kappa_zz(E[ii], ISM_phase[ni], mass = cst.mp, kmin = (50*cst.pc)**(-1), q = 5./3, I = 1e-2)
    
    if (E[ii] < max(emax_SNR)) : 
        for xi in range(1, len(SNR.get("R_SNR"))) : 
            if (SNR.get("t_SNR")[xi] > tesc_SNR[ii] >= SNR.get("t_SNR")[xi-1]) : 
                Resc[ii] = 0.5*(SNR.get("R_SNR")[xi] + SNR.get("R_SNR")[xi-1])


# {"R_SNR"   :r_new, 
#             "t_SNR"   :t_new, 
#             "u_sh"    :u_sh, 
#             "t_SED"   :tfree, 
#             "t_PDS"   :tPDS,
#             "t_MCS"   :tMCS,
#             "t_merge" :tmerge,
#             "t_max"   : tmax, 
#             "r_SED"   : R_free,
#             "r_PDS"   : R_PDS,
#             "r_MCS"   : R_MCS,
#             "r_merge" : R_merge}




size_x = 4
size_y = 3
sub_x  = 2
sub_y  = 2
fig = plt.figure(figsize=(size_x*sub_x,size_y*sub_y))

gs = gridspec.GridSpec(ncols= sub_x, nrows = sub_y, figure = fig )
gs.update(wspace=0.5, hspace=0.5) # set the spacing between axes.


ax00l = fig.add_subplot(gs[0, 0])
ax00r = ax00l.twinx()
ax00l.loglog(SNR.get("t_SNR")/cst.kyr, SNR.get("R_SNR")/cst.pc, lw = 2, c = "blue")
ax00r.loglog(SNR.get("t_SNR")/cst.kyr, SNR.get("u_sh")/cst.kms, lw = 2, c = "red")

ax00l.scatter(SNR.get("t_SED")/cst.kyr, SNR.get("r_SED")/cst.pc, marker="o",c="black",label="$t_\\mathrm{Sed}$")
ax00l.scatter(SNR.get("t_PDS")/cst.kyr, SNR.get("r_PDS")/cst.pc, marker="s",c="black",label="$t_\\mathrm{PDS}$")
# ax00l.plot(SNR.get("t_MCS")/cst.kyr, SNR.get("r_MCS")/cst.pc, marker="v",c="black",label="$t_\\mathrm{MCS}$")

ax00l.set_ylabel("$R_\\mathrm{SNR}$ [pc]", color = "blue")
ax00r.set_ylabel("$v_\\mathrm{sh}$ [km/s]", color = "red")
ax00l.set_xlabel("$t$ [kyr]")
ax00l.set_xlim(1e-1, 1e3)
ax00l.set_ylim(1e0, 1e2)
ax00l.legend(loc = "best", ncol = 2)

ax01l = fig.add_subplot(gs[0, 1]) 
ax01r = ax01l.twinx() 
ax01l.loglog(E/cst.GeV, tesc_SNR/cst.kyr, lw = 2, c = "red")
ax01r.semilogx(E/cst.GeV, Resc/cst.pc, lw = 2, c = "blue")

ax01r.set_ylabel("$R_\mathrm{esc}$ [pc]", color = "blue") 
ax01l.set_ylabel("$t_\mathrm{esc}$ [kyr]", color = "red") 
ax01l.set_xlabel("$E$ [GeV]")
ax01r.set_ylim(3, 15)
ax01l.set_xlim(1e-1, 1e5)





ax11 = fig.add_subplot(gs[1, :]) 
for ni in range(len(ISM_phase)) : 
    ax11.loglog(E/cst.GeV, kzz_d[ni], c = ISM_color[ni], ls = '-', lw = 2, label = ISM_name[ni])
ax11.loglog(E/cst.GeV, kzz_i, c = "black", ls = "--", lw = 2, label = "ISM independant")

ax11.legend(ncol = len(ISM_name)+1, loc = "best")

ax11.set_ylabel("$D_0$ [cm$^2$/s]")
ax11.set_xlabel("$E$ [GeV]")
ax11.set_xlim(1e-1, 1e5)


plt.savefig("all_in_one.pdf", bbox_inches="tight")



























