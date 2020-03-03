#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 13:55:50 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import mathmethods as math
#import phases_collection as ism
import constants as cst





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







# Ion-neutral damping of Alfven waves 
# This routine return the solution of the Alfven waves 
# dispersion relation and the Alfven velocity
# Variables : 
# - position_index : refers to the xi index position in the grid
# - E refers to the reseonnant CRs energy in erg
# - medium_props : correspond to a dictionnary containing arrays of ism props as energy and position
def indamping_alfven(position_index , E, medium_props) : 
    # Import values of phases_collection
    T  = medium_props.get('T')[position_index] 
    B  = medium_props.get('B')[position_index] 
    mn = medium_props.get('mn')[position_index] 
    mi = medium_props.get('mi')[position_index] 
    X  = medium_props.get('X')[position_index] 
    ni = medium_props.get('ni')[position_index] 
    nn = medium_props.get('nn')[position_index] 
    
#    print "T = ",T," K"
#    print "B = ",B," G"
#    print "mn = ",mn/cst.mp," mp"
#    print "mi = ",mi/cst.mp," mp"
#    print "X = ",X
#    print "ni = ",ni," cm^-3"
#    print "nn = ",nn," cm^-3"
    
    # Some relations 
    rl = (E)/(cst.e*B)
    k = 1/rl
#    chi = (mn/mi)*(X**(-1.)-1.)
    rho_n = mn*nn
    rho_i = mi*ni
    chi = (mn*nn)/(mi*ni)
    VAi = B/np.sqrt(4*np.pi*mi*ni)
    if (T <= 50) : 
#        nu_in = nn*2.1e-9
        nu_in = 2*nn*8.4e-9*(50/1e4)**0.4
    if (T > 50) : 
        nu_in = 2*nn*8.4e-9*(T/1e4)**0.4
    nu_ni = chi**(-1.)*nu_in
    
#    print "E = ",E/cst.TeV," TeV"
#    print "rl = ",rl/cst.pc," pc"
#    
#    print "chi = ",chi
#    print "Vai = ",VAi," cm/s"
#    print "nu_in = ",nu_in," s^-1"
#    print "nu_ni = ",nu_ni," s^-1"
    # We set these values to zero for the moment
    nu_n  = 0.
    theta = 0.
    # We calculate these values 
    a =  k**2*nu_n + (1 + chi)*nu_ni
    b = (1/4.)*(k**2*np.cos(theta)**2*VAi**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1 + chi)*nu_ni)**2)
    c = (1/8.)*((k**2*nu_n + (1 + chi)*nu_ni)*chi*k**2*nu_n*nu_ni + chi*nu_ni*k**2*np.cos(theta)**2*VAi**2)
    wi = math.cardano3(a, b, c)
    wr = np.sqrt(3*wi**2 + 2*(k**2*nu_n + (1 + chi)*nu_ni)*wi + k**2*np.cos(theta)**2*VAi**2 + chi*k**2*nu_n*nu_ni)
    
    if (np.isnan(wr)) : 
        wr = 0.
    
    return dict(wr=wr, wi=-wi, VA=wr/k, nu_ni=nu_ni, chi=chi, rho_n=rho_n, rho_i=rho_i)






def indamping_alfven_nopos(E, medium_props) : 
    # Import values of phases_collection
    T  = medium_props.get('T')
    B  = medium_props.get('B')
    mn = medium_props.get('mn')
    mi = medium_props.get('mi') 
    X  = medium_props.get('X')
    ni = medium_props.get('ni')
    nn = medium_props.get('nn')
    
#    print "T = ",T," K"
#    print "B = ",B," G"
#    print "mn = ",mn/cst.mp," mp"
#    print "mi = ",mi/cst.mp," mp"
#    print "X = ",X
#    print "ni = ",ni," cm^-3"
#    print "nn = ",nn," cm^-3"
    
    # Some relations 
    rl = (E)/(cst.e*B)
    k = 1/rl
    chi = (mn/mi)*(X**(-1.)-1.)
#    rho_n = mn*nn
#    rho_i = mi*ni
#    chi = (mn*nn)/(mi*ni)
    VAi = B/np.sqrt(4*np.pi*mi*ni)
    if (T <= 50) : 
#        nu_in = nn*2.1e-9
        nu_in = 2*nn*8.4e-9*(50/1e4)**0.4
    if (T > 50) : 
        nu_in = 2*nn*8.4e-9*(T/1e4)**0.4
    nu_ni = chi**(-1.)*nu_in
    
#    print "E = ",E/cst.TeV," TeV"
#    print "rl = ",rl/cst.pc," pc"
#    
#    print "chi = ",chi
#    print "Vai = ",VAi," cm/s"
#    print "nu_in = ",nu_in," s^-1"
#    print "nu_ni = ",nu_ni," s^-1"
    # We set these values to zero for the moment
    nu_n  = 0.
    theta = 0.
    # We calculate these values 
    a =  k**2*nu_n + (1 + chi)*nu_ni
    b = (1/4.)*(k**2*np.cos(theta)**2*VAi**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1 + chi)*nu_ni)**2)
    c = (1/8.)*((k**2*nu_n + (1 + chi)*nu_ni)*chi*k**2*nu_n*nu_ni + chi*nu_ni*k**2*np.cos(theta)**2*VAi**2)
    wi = math.cardano3(a, b, c)
    wr = np.sqrt(3*wi**2 + 2*(k**2*nu_n + (1 + chi)*nu_ni)*wi + k**2*np.cos(theta)**2*VAi**2 + chi*k**2*nu_n*nu_ni)
    
    if (np.isnan(wr)) : 
        wr = 0.
    
    return dict(wr=wr, wi=-wi, VA=wr/k, nu_ni=nu_ni, chi=chi)


# Non-correlated interactions with large scale turbulence 

def damping_lazarian(position_index, E, medium_props) : 
    # Import values of phases_collection
    T  = medium_props.get('T')[position_index] 
    B  = medium_props.get('B')[position_index] 
    mn = medium_props.get('mn')[position_index] 
    mi = medium_props.get('mi')[position_index] 
    X  = medium_props.get('X')[position_index] 
    ni = medium_props.get('ni')[position_index] 
    nn = medium_props.get('nn')[position_index] 
    
    alfven_infos = indamping_alfven(position_index , E, medium_props)
    ca = alfven_infos.get('VA')
    
    if (ca < 1e-30) : 
        damp_lz = 0.
    
    L = 50*cst.pc # Turbulence scale injection
    rl = (E)/(cst.e*B)
#    rl = E*1e9/B/300
    

    # Some relations 
    eps   = ni*mi/(nn*mn)
    cai   = B/np.sqrt(4*np.pi*mi*ni)
    can   = B/np.sqrt(4*np.pi*(mi*ni + mn*nn))
    if (T > 100) : coll = nn*2*8.4e-9*(T/1e4)**(0.4)
    if (T<= 100) : coll = nn*2.1e-9
    gamma_ad = 5./3
    mu_n  = mn/cst.mp
    nu_nn = 1.5e-10*T**(0.31)
    c_n   = 9.79e5*np.sqrt(gamma_ad/mu_n)*(T*cst.kbolz*6.242e11)**(0.5)    
    xi_n  = (nn*mi)/(nn*mn + ni*mi)
    nu_ndamp = 0
    ln = c_n*(eps/coll + 1./(nn*nu_nn))
    if (X < 0.2) : nu_n = nn/(nn + ni)*c_n*ln
    if (X >= 0.2) : nu_n = 0
    
    VL = max(can, c_n)
    MA = VL/can
    lA = L*MA**(-3.)
    Lcm = L
    
    if (nu_n == 0) : 
        return 0 
    
    kdampar_sub = (-(nu_n + pow(cai,2.)/coll) + pow(pow(nu_n + pow(cai,2.)/coll,2) + 8*can*nu_n*Lcm*pow(MA,-4.)/xi_n, 0.5))/(2*nu_n*Lcm*pow(MA,-4.))
    lmin_sub = pow(kdampar_sub*np.sqrt(1 + Lcm*pow(MA,-4.)*kdampar_sub),-1.)
    
    kdampar_super = (-(nu_n + pow(cai,2.)/coll) + pow(pow(nu_n + pow(cai,2.)/coll,2) + 8*can*nu_n*lA/xi_n, 0.5))/(2*nu_n*lA)
    lmin_super = pow(kdampar_super*np.sqrt(1 + lA*kdampar_super),-1.)
    
    
    
    if (MA <= 1.) : 
        if (rl < pow(lmin_sub,4./3)*pow(MA,4./3)*pow(Lcm,-1./3)) : 
            damp_lz = 0.
        if (rl >= pow(lmin_sub,4./3)*pow(MA,4./3)*pow(Lcm,-1./3) and rl < Lcm*pow(MA,4)) : 
            damp_lz = -ca*pow(MA,2.)*pow(rl, -1./2)*pow(Lcm,-1./2)
        if (rl >= Lcm*pow(MA,4.) and rl < Lcm*MA) : 
            damp_lz = -ca*pow(MA, 8./3)*pow(rl, -2./3)*pow(Lcm, -1./3)
        if (rl >= Lcm*MA) : 
            if (rl < Lcm) : 
                damp_lz = -pow(MA,2.)*ca/Lcm
            if (rl >= Lcm) : 
                damp_lz = -pow(MA,2)*ca*Lcm*pow(rl,-2.)
    if (MA > 1.) : 
        if (rl < pow(lmin_super,4./3)*MA*pow(Lcm,-1./3)) : 
            damp_lz = 0.
        if (rl >= pow(lmin_super,4./3)*MA*pow(Lcm,-1./3) and rl < Lcm*pow(MA,-3.)) :
            damp_lz = -ca*pow(MA,3./2)*pow(Lcm,-1./2)*pow(rl,-1./2)
        if (rl >= Lcm*pow(MA,-3.)) : 
            damp_lz = -ca*MA*pow(Lcm,-1./3)*pow(rl,-2./3)
    
    return damp_lz



def non_linear_landau_damping(T, Ip, Im, mi, q, B0, Ecr) : 
    vi = np.sqrt(cst.kbolz*T/mi)
    dB2 = Ip**2 + Im**2
    kr  = q*B0/Ecr
    gamma_nlld = - np.sqrt(np.pi/8.)*vi*kr*dB2
    return gamma_nlld
