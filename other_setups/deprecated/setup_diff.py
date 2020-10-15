#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:17:05 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import sys 
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from scipy import interpolate


# We refers to the tools setup folder
# import sys
import os 
import shutil 
sys.path.append('../tools/')
 
import freader as fr 
import fwritter as fw
import constants as cst
import d1_grid_generator as grid 
import phases_collection as ism
import damping as dp
import background_diffusion_coefficient as bdc



###############################################################################
#   INFORMATIONS COMING FROM THE NAMELIST                                     #
###############################################################################
# !!! Do not modify this part !!! #
import namelist_diff as nml

# Phase space 
X         = nml.X
E         = nml.E
nx        = nml.NX
ne        = nml.NE

x_center = nml.box_center
x_center_index = int(x_center/(X[1] - X[0]))


# ISM secondary variables 
ism_values = nml.ism_values
B  = ism_values.get("B")
nn = ism_values.get("nn")
ni = ism_values.get("ni")
mn = ism_values.get("mn")
mi = ism_values.get("mi")
T  = ism_values.get("T")
Xi = ism_values.get("X")
va = ism_values.get("VA")
###############################################################################



###############################################################################
#   INITIAL ISM CONDITIONS                                                    #
###############################################################################
# You can modify this part #

# Initial diffusion coefficients -> Initial waves I
d00      = np.zeros(len(X)) # Spatial component of the background diffusion coefficient of CRs 
D        = np.zeros((len(E), len(X))) # Initial diffusion coefficient of CRs 
Db       = np.zeros((len(E), len(X))) # Bohm diffusion coefficient of CRs 
Ip       = np.zeros((len(E), len(X))) # Initial resonnant fowarding waves energy density 
Im       = np.zeros((len(E), len(X))) # Initial resonnant backwarding waves energy density 

# Calculation of the Alfven velocity and the ion-neutral damping 
VA             = np.zeros((len(E), len(X)))
gamma_in       = np.zeros((len(E), len(X)))
gamma_lazarian = np.zeros((len(E), len(X)))
gamma_nlld     = np.zeros((len(E), len(X)))
gamma_tot      = np.zeros((len(E), len(X)))
Pcr            = np.zeros((len(E), len(X)))
Pe             = np.zeros((len(E), len(X))) # Electrons pressure distribution 



def door(X, X1, X2, V) : 
    if (X > X1 and X < X2) : 
        return V
    else : 
        return 0. 



for ei in range(len(E)) : 
    print (ei,"over ",len(E))
    for xi in range(len(X)) : 
        # VA[ei][xi] = 1.*cst.pc/cst.yr
        if (X[xi] >= x_center) : 
            VA[ei][xi] = 0.5*cst.pc/cst.yr
        else : 
            VA[ei][xi] = -0.5*cst.pc/cst.yr
        
        Pcr[ei][xi] = 0. + door(X[xi], 950.*cst.pc, 1050.*cst.pc, 1.)
        
        D[ei] = (5*cst.pc)**2/cst.yr
        Db[ei][xi] = 1.
        Ip[ei][xi] = 10**(-29)
        Im[ei][xi] = 0.



# for xi in range(len(X)) : 
#     # Normalisation of the background diffusion coefficient
#     d00[xi] = 1e28

# for e in range(len(E)) : 
#     if (e % 1 == 0) : 
#         print ("Background diffusion model : "+str(round(e/len(E)*100.,2))+" %")
#     for xi in range(len(X)) : 
#         # Bohm diffusion coefficient
#         Db[e][xi] = (4*cst.c)/(3*np.pi)*(E[e]/(cst.e*B[x_center_index])) 
        
#         # Background diffusion coefficient
#         if (nml.bdiff_model == "ISM_independant") : 
#             D[e][xi]  = d00[xi]*(E[e]/(10.*cst.GeV))**0.5                      # ISM independant method
#         if (nml.bdiff_model == "ISM_dependant") : 
#             medium_props = {"B"  : B[xi], 
#                             "mn" : mn[xi],
#                             "nn" : nn[xi],
#                             "mi" : mi[xi],
#                             "ni" : ni[xi],
#                             "X"  : Xi[xi],
#                             "T"  : T[xi]}
#             D[e][xi]  = bdc.Kappa_zz(E[e], 
#                                      medium_props, 
#                                      mass = cst.mp, 
#                                      kmin = (50.*cst.pc)**(-1), 
#                                      q = 5./3, 
#                                      I = 1e-4)                              # ISM dependant method
#         # Background rates of turbulence
#         Ip[e][xi] = Db[e][xi]/D[e][xi]*0.5
#         Im[e][xi] = Db[e][xi]/D[e][xi]*0.5



# for e in range(len(E)) : 
#     if (e % 10 == 0) : 
#         print ("Waves damping and initial CRs distributions : "+str(round(e/len(E)*100.,2))+" %")
#     for xi in range(len(X)) : 
#         in_damping = dp.indamping_alfven(xi , E[e], ism_values) 
        
#         # VA is a vectorial field
#         # Advection append from the center of the source to the edges of the simulation
#         if (X[xi] >= x_center) : 
#             VA[e][xi] = +va[e][xi]
#         else : 
#             VA[e][xi] = -va[e][xi]

#         if (nml.in_damping) : 
#             gamma_in[e][xi]       = in_damping.get('wi')
#         if (nml.lz_damping) : 
#             gamma_lazarian[e][xi] = -dp.damping_lazarian(xi, E[e], ism_values)
#         if (nml.nlld_damping) : 
#             gamma_nlld[e][xi]     = -dp.non_linear_landau_damping( T[xi], Ip[e][xi], Im[e][xi], 
#                                                                   mi[xi],     cst.e, 
#                                                                    B[xi],      E[e])
#         gamma_tot[e][xi]          = gamma_in[e][xi] + gamma_lazarian[e][xi] + gamma_nlld[e][xi]
#         # Background CRs spectra
#         Pcr[e][xi]                = nml.Pcr_1GeV*(E[e]/cst.GeV)**(-2.7)
#         Pe[e][xi]                 = nml.Pe_1GeV*(E[e]/cst.GeV)**(-3.1)








###############################################################################
#    END : INITIAL ISM CONDITIONS                                             #
###############################################################################

###############################################################################
#      WRITE THE INITAL CONDITIONS                                            #
###############################################################################
print ("Writing Initial Conditions")
# !!! Do not modify this part !!! #
fw.write1D( "Alfven.dat", nx=nx, ne=ne, variable=VA.T, path=nml.total_path+"/data_ini/") 
fw.write1D("DCRpara.dat", nx=nx, ne=ne,  variable=D.T, path=nml.total_path+"/data_ini/")
fw.write1D(  "DBohm.dat", nx=nx, ne=ne, variable=Db.T, path=nml.total_path+"/data_ini/")
fw.write1D(     "Ip.dat", nx=nx, ne=ne,  variable=Ip.T, path=nml.total_path+"/data_ini/")   
fw.write1D(     "Im.dat", nx=nx, ne=ne,  variable=Im.T, path=nml.total_path+"/data_ini/")  
fw.write1D(    "Pcr.dat", nx=nx, ne=ne,variable=Pcr.T, path=nml.total_path+"/data_ini/")    
fw.write1D(    "Pe.dat", nx=nx, ne=ne,variable=Pe.T, path=nml.total_path+"/data_ini/")   
fw.write1D("damping.dat", nx=nx, ne=ne, variable=gamma_tot.T, path=nml.total_path+"/data_ini/")  
fw.write1Daxis("X.dat", variable=X, nx=nx, path=nml.total_path+"/data_ini/")
fw.write1Daxis("E.dat", variable=E, nx=ne, path=nml.total_path+"/data_ini/")
fw.write1Daxis("B.dat", variable=B, nx=nx, path=nml.total_path+"/data_ini/")

variables = {"NX"       : nx,
             "NE"       : ne,
             "ni"       : ni[x_center_index],
             "X"        : Xi[x_center_index],
             "mn"       : mn[x_center_index],
             "T"        : T[x_center_index],
             "center"   :x_center,
             "center_index":x_center_index,
             "B"        : B[x_center_index]}

fw.fileWrite("parameters", variables = variables, path=nml.total_path+"/", ext='.dat') 
###############################################################################


###############################################################################
#      COPY TOOLS AND ANALYSIS FILES                                          #
###############################################################################
# !!! Do not modify this part !!! #
shutil.copy("../data_analysis/show_data.py", nml.total_path+"/data_analysis/") 
shutil.copy("../data_analysis/pcr_ip_2D.py", nml.total_path+"/data_analysis/") 

shutil.copy("../tools/freader.py", nml.total_path+"/tools/") 
shutil.copy("../tools/constants.py", nml.total_path+"/tools/") 










  