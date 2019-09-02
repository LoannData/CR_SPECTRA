#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 08:26:46 2019

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import mathmethods as math
import phases_collection as ism
import constants as cst

from scipy.interpolate import interp1d

def collision_rate(specie1, specie2, phase) : 
    if (specie1 == "ion" and specie2 == "neutral"): 
        if (phase.get("T")<= 50) : 
            nu_in = phase.get("nn")*2.1e-9
        else : 
            nu_in = 2*phase.get('nn')*8.4e-9*(phase.get("T")/1e4)**0.4
        return nu_in
    if (specie1 == "neutral" and specie2 == "ion") : 
        chi = (phase.get("mn")/phase.get("mi"))*(phase.get("X")**(-1.)-1.)
        nu_in = collision_rate("ion", "neutral", phase)
        nu_ni = chi**(-1.)*nu_in
        return nu_ni 

def cr_escape_time_model(option, model, Ecr, props) : 
    if (option == 2 and model == "model1") : 
        c1 = props.get("c1")
        c2 = props.get("c2")
        return c2*(Ecr/c1)**(-0.5)
    if (option == 2 and model == "N16") : 
        phase = props.get("phase")
        spec_index = props.get("cr_index")
        E = np.array([10.*cst.GeV, 30.*cst.GeV, 100.*cst.GeV, 300.*cst.GeV, 1e3*cst.GeV, 3e3*cst.GeV, 10e3*cst.GeV, 30e3*cst.GeV])
        if (phase == "WNM") : 
            if (spec_index == 2.2) : 
                t = np.array([38.68850, 16.354, 6.838483, 3.2713, 1.4532, 0.704590, np.nan, 0.207302])*cst.kyr 
            if (spec_index == 2.0) : 
                t = np.array([26.529, 12.620, 6.02296, 3.19788, 1.64029, 0.8756588,     np.nan, 0.330993])*cst.kyr 
        if (phase == "CNM") : 
            if (spec_index == 2.2) : 
                t = np.array([8.24146, 3.652, 1.22933, np.nan, 1.169809, 0.718255, 0.1, np.nan])*cst.kyr
            if (spec_index == 2.0) : 
                t = np.array([5.4321, 2.6174, 1.0571,    np.nan, 1.34488,  1.3592, 0.25691,    0.043])*cst.kyr
        if (phase == "DiM") : 
            if (spec_index == 2.2) : 
                t = np.array([2.964, np.nan, np.nan, 3.11714, 2.8569, 0.53177, 0.07047, 0.02210])*cst.kyr
            if (spec_index == 2.0) : 
                t = np.array([1.8941,    np.nan,    np.nan, 3.0516, 3.3622,  1.0786, 0.17705, 0.04029])*cst.kyr
        Elin = np.zeros(len(E))
        tlin = np.zeros(len(t))
        for ii in range(len(t)) : 
            Elin[ii] = np.log10(E[ii])
            tlin[ii] = np.log10(t[ii])
        idx = np.isfinite(Elin) & np.isfinite(tlin)
        
        p = interp1d(Elin[idx], tlin[idx], fill_value="extrapolate", kind="cubic")
#        z = np.polyfit(Elin[idx], tlin[idx], 4)
#        p = np.poly1d(z)
        if (Ecr/cst.GeV < 5. or Ecr/cst.TeV > 30.) : print "Take care ! The escape model N16 is only available between 5 GeV and 30 TeV"
        return 10**(p(np.log10(Ecr)))
        

def cr_escape_radius_model(option, model, time, props) : 
    if (option == 2 and model == "snr1") : 
        return 25*cst.pc*(time/(25*cst.kyr))**(0.5)

#    xi_n = 1 # = 1 For solar abundances
#    phi_c = 1 # Actual thermal cond. / the Sptitzer (1962) value 
#    beta = 2. 
#    vej8 = 10*(Esn/Mej)**(1./2)
#    C06 = 9.79e5*np.sqrt(5./3/(ph_mn[index]/mp))*np.sqrt(ph_T[index]*kbolz*6.24e11)/1e6 
#    tfree = 0.3*Esn**(-1./2)*Mej*n**(-1./3) # [kyr]
#    tsf = 3.61e4*(Esn**(3./14))/(n**(4./7))*1e-3 # [kyr]
#    tPDS = tsf/np.exp(1) # [kyr]
#    tMCS = min(61*vej8**3/(xi_n**(9./14)*n**(3./7)*Esn**(3./14)), 476/(xi_n*phi_c)**(9./14))*tPDS # [kyr]
#    tmerge = 153*(Esn**(1./14)*n**(1./7)*xi_n**(3./14)/(beta*C06))**(10./7)*tPDS # [kyr]
#    
#    ctimes = [tfree*1e3, tPDS*1e3, tMCS*1e3, tmerge*1e3]
#    
#    
#    
#    RMCS = (4.66*(tMCS)*(1 - 0.939*(tMCS)**(-0.17) + 0.153*(tMCS)**(-1.)))**(1./4)
    if (option == 2 and model == "C88&TM99") : 
        phase = props.get("phase")
        ni = phase.get("ni")
        fion = phase.get("X")
        mn = phase.get("mn")
        T = phase.get("T")
        Esn = props.get("Energy_SNR")
        Mej = props.get("Mass_ejecta")
        Esn = Esn/1e51
        Esn51 = Esn
        Mejs  = Mej 
        n = ni/fion  # Total density of the medium [cm^-3]
        xi_n = 1 # = 1 For solar abundances
        phi_c = 1 # Actual thermal cond. / the Sptitzer (1962) value 
        beta = 2. 
        vej8 = 10*(Esn/Mej)**(1./2)
        C06 = 9.79e5*np.sqrt(5./3/(mn/cst.mp))*np.sqrt(T*cst.kbolz*6.24e11)/1e6 
        tfree = 0.3*Esn**(-1./2)*Mej*n**(-1./3) # [kyr]
        tsf = 3.61e4*(Esn**(3./14))/(n**(4./7))*1e-3 # [kyr]
        tPDS = tsf/np.exp(1) # [kyr]
        tMCS = min(61*vej8**3/(xi_n**(9./14)*n**(3./7)*Esn**(3./14)), 476/(xi_n*phi_c)**(9./14))*tPDS # [kyr]
        tmerge = 153*(Esn**(1./14)*n**(1./7)*xi_n**(3./14)/(beta*C06))**(10./7)*tPDS # [kyr]
        def RSNR_ST(timekyr) : 
            return  5.0*(Esn51/n)**(1./5)*(1 - (0.05*Mejs**(5./6)/(Esn51**(1./2)*n**(1./3)*timekyr)))*(2./5)*(timekyr)**(2./5)
        
        def RSNR_PDS(timekyr) : 
            R0 = RSNR_ST(tPDS)
            return R0*(timekyr/tPDS)**(3./10)
        
        def RSNR_MCS(timekyr) : 
            R0 = RSNR_PDS(tMCS)
            return R0*(timekyr/tMCS)**(1./4)
        
        def RSNR_free(timekyr) : 
            R0 = RSNR_ST(tfree)
            return R0*(timekyr/tfree)**(1.)
        # RNSR is in [pc]
        if (time/cst.kyr < tfree) : 
            RSNR = RSNR_free(time/cst.kyr)
        if (time/cst.kyr > tfree and time/cst.kyr < tPDS) : 
            RSNR = RSNR_ST(time/cst.kyr)
        if (time/cst.kyr > tPDS and time/cst.kyr < tMCS) : 
            RSNR = RSNR_PDS(time/cst.kyr)
        if (time/cst.kyr > tMCS and time/cst.kyr < tmerge) : 
            RSNR = RSNR_MCS(time/cst.kyr)
        if (time/cst.kyr > tmerge) : 
            RSNR = np.NaN
#        print "###########################"
#        print "time = ",time/cst.kyr," kyr"
#        print "tfree = ",tfree," kyr"
#        print "tPDS = ",tPDS," kyr"
#        print "tMCS = ",tMCS," kyr"
#        print "tmerge = ",tmerge," kyr"
        return 4.*RSNR*cst.pc 
        

        

        





# Test for cr_escape_time_model N16
#####################################################
#option = 2 
#model = "N16"
#props = {"phase":"CNM", "cr_index":2.2}
#Ecr = np.logspace(np.log10(5.*cst.GeV), np.log10(50.*cst.TeV), 100)
#tesc = np.empty(len(Ecr))
#
#for ii in range(len(Ecr)) : 
#    tesc[ii]= cr_escape_time_model(option, model, Ecr[ii], props)
#
#plt.loglog(Ecr/cst.GeV, tesc/cst.kyr, c="red")
#####################################################

# Test for SNR evolution model C88&TM99
#####################################################
#option = 2 
#model = "C88&TM99"
#props_Resc = {"phase":ism.WNM, "Energy_SNR":1e51, "Mass_ejecta":1.}
#time = np.logspace(np.log10(0.01*cst.kyr), np.log10(1000*cst.kyr), 100)
#RSNR = np.empty(len(time))
#
#for ii in range(len(time)) : 
#    RSNR[ii] = cr_escape_radius_model(option, model, time[ii], props_Resc)
#
#plt.loglog(time/cst.kyr, RSNR/cst.pc)
#####################################################