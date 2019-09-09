#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:27:43 2019

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 

import sys 
sys.path.append('../tools/')
import freader as fr 
import constants as cst


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

def PlotXSpace(var, file_number, energy, fig_size_x=10, fig_size_y=6, xlabel="z [pc]", ylabel="None",
               savefig=False, figname="None") : 
    X, E, data = getData(var, file_number)
    
    plt.figure(figsize=(fig_size_x, fig_size_y))
    
    
    for ei in range(len(energy)) : 
        e_id = 0
        closest = np.inf 
        for ii in range(len(E)) : 
            if (abs(E[ii]-energy[ei]) < closest) : 
                e_id = ii
                closest = abs(E[ii]-energy[ei])
    
        plt.loglog(np.array(X)/cst.pc, data.T[e_id], 
                   label="$E_\\mathrm{CR} = $"+str(round(E[e_id]/cst.GeV,2))+" GeV")
    
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    
    if (savefig) : 
        plt.savefig("./"+figname+".pdf")

def PlotESpace(var, file_number, position, fig_size_x=10, fig_size_y=6, xlabel="E [GeV]", ylabel="None",
               savefig=False, figname="None") : 
    X, E, data = getData(var, file_number)
    
    plt.figure(figsize=(fig_size_x, fig_size_y))
    
    
    for ei in range(len(position)) : 
        e_id = 0
        closest = np.inf 
        for ii in range(len(X)) : 
            if (abs(X[ii]-position[ei]) < closest) : 
                e_id = ii
                closest = abs(X[ii]-position[ei])
    
        plt.loglog(np.array(E)/cst.GeV, data[e_id], 
                   label="$X = $"+str(round(X[e_id]/cst.pc,2))+" pc")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    
    if (savefig) : 
        plt.savefig("./"+figname+".pdf")



###############################################################################
# Put your instructions here !                                                #
###############################################################################
PlotXSpace("Dcr", 0, np.logspace(np.log10(100.*cst.GeV), np.log10(10.*cst.TeV), 4),
           ylabel="$P_\\mathrm{CR}$ [erg/cm$^3$]")

PlotESpace("Dcr", 0, np.logspace(np.log10(50.*cst.pc), np.log10(500.*cst.pc), 4),
           ylabel="$P_\\mathrm{CR}$ [erg/cm$^3$]")






#plt.figure(figsize=(6, 8))
##plt.loglog(E, data[500])
#plt.loglog(np.array(X)/3.086e18, data.T[1])
#plt.show()