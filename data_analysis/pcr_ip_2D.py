#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 09:25:14 2019

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

import sys 
sys.path.append('../tools/')
import freader as fr 
import constants as cst

# Functions to read the data
###############################################################################
def readDataXE(file_name, NX, NE) : 
    e_cols = []
    for ei in range(NE) : 
        e_cols.append(str(ei))
    df = pd.read_csv(file_name, usecols=e_cols)
    var = np.empty((NX, NE))
    for xi in range(NX) : 
        for ei in range(NE) : 
            var[xi][ei] = df[str(ei)].iloc[xi]
    return var

def readAxis(file_name) : 
    myfile = open(file_name,"r").readlines()
    for ii in range(len(myfile)) : 
        myfile[ii] = myfile[ii].replace("\n","")
        myfile[ii] = float(myfile[ii])
    return myfile


def getData(var, file_number) : 
    if (file_number < 10) : 
        sfile_number = "000"+str(file_number)
    if (file_number >= 10 and file_number < 100) : 
        sfile_number = "00"+str(file_number)
    if (file_number >= 100 and file_number < 1000) : 
        sfile_number = "0"+str(file_number)
    if (file_number >= 1000) : 
        sfile_number = str(file_number)
    nx = fr.search("../parameters.dat", "NX")
    ne = fr.search("../parameters.dat", "NE")
    X = readAxis("../data_ini/X.dat")
    E = readAxis("../data_ini/E.dat")
    if (var == "Pcr" or var == "Ip") : 
        data = readDataXE("../data_out/"+var+"_"+sfile_number+".dat", 2**int(nx), 2**int(ne))
    if (var == "Dcr") : 
        data = np.empty((len(X), len(E)))
        Ip = readDataXE("../data_out/"+"Ip"+"_"+sfile_number+".dat", 2**int(nx), 2**int(ne))
        Db = readDataXE("../data_ini/"+"DBohm"+".dat", 2**int(nx), 2**int(ne))
        for ii in range(len(Ip)) : 
            for jj in range(len(Ip[ii])) : 
                data[ii][jj] = Db[ii][jj]/Ip[ii][jj]
    return X, E, data

def getTimeID(time_test, delta_t, t_ini, kind="linear") : 
    if (kind == "linear") : 
        out_id = int(round((time_test - t_ini)/delta_t , 0)) 
    return out_id 
###############################################################################
    

t_ini   = 0.*cst.kyr  # Initial output plot time 
t_max   = 72.*cst.kyr # Max output time of the simulation 
delta_t = 1.*cst.kyr  # Time distance between two plots 



time_test = 35.*cst.kyr # We choose a fixed output time 
out_id = getTimeID(time_test, delta_t, t_ini) # We convert it in an index number 

z_wnmcnm = 100.*cst.pc
z_wnmcnm = np.log10(z_wnmcnm/cst.pc)

# We get our data 
X, E, Pcr = getData("Pcr", out_id)
X, E, Ip  = getData("Dcr",  out_id)




EV, XV = np.meshgrid(E, X, sparse=False, indexing='xy')
cmap = plt.get_cmap('PiYG')


XV_log = np.log10(XV/cst.pc)
EV_log = np.log10(EV/cst.GeV)

XV_log[0] = -100

fig, (ax0, ax1) = plt.subplots(ncols = 2, sharey=False, figsize=(16, 6))
#fig.figure(figsize=(6, 10))

im0 = ax0.pcolormesh(XV_log, EV_log, np.log10(Pcr), cmap=cmap)
ax0.set_ylim(np.log10(E[0]/cst.GeV), np.log10(E[-1]/cst.GeV))
ax0.set_xlim(np.log10(X[1]/cst.pc), np.log10(X[-1]/cst.pc))
ax0.axvline(z_wnmcnm, c="black")
fig.colorbar(im0, ax=ax0, cax = fig.add_axes([0.05, 0.1, 0.02, 0.8]), 
             label = "$\\log(P_\\mathrm{CR})$ [erg/cm$^3$]")

im1 = ax1.pcolormesh(XV_log, EV_log, np.log10(Ip), cmap=cmap)
ax1.set_xlim(np.log10(X[1]/cst.pc), np.log10(X[-1]/cst.pc))
ax1.set_ylim(np.log10(E[0]/cst.GeV), np.log10(E[-1]/cst.GeV))
ax1.axvline(z_wnmcnm, c="black")
#cb1 = fig.colorbar(im1, ax=ax1)
fig.colorbar(im1, ax=ax1, cax = fig.add_axes([0.92, 0.1, 0.02, 0.8]), 
             label = "$\\log(D_\\mathrm{CR})$ [cm$^2$/s]")


ax0.set_ylabel("E [GeV]")
ax0.yaxis.set_label_position("right")
ax0.yaxis.tick_right()
ax1.yaxis.tick_left()


ax0.set_xlabel("$\\log(z)$ [pc]")
ax1.set_xlabel("$\\log(z)$ [pc]")

fig.suptitle("Time = "+str(round(time_test/cst.kyr,2))+" kyrs", fontsize=16)


#plt.pcolormesh(XV_log, EV_log, np.log10(Pcr), cmap=cmap)
#plt.ylim(np.log10(E[0]/cst.GeV), np.log10(E[-1]/cst.GeV))
#plt.xlim(np.log10(X[1]/cst.pc), np.log10(X[-1]/cst.pc))











