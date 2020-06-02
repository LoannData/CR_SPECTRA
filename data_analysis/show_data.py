#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:24:07 2020

@author: lbrahimi
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
import pandas as pd 

import sys 
sys.path.append('../tools/')
import freader as fr 
import constants as cst

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def colorArray(c1, c2, n = 10) : 
    colors = []
    for ii in range(n) : 
        colors.append(colorFader(c1, c2, mix = ii/(n-1)))
    return colors 

###############################################################################
# DATA READING ROUTINES                                                       #
###############################################################################
# !!! DO NOT MODIFY !!!                                                       # 

def readDataXE(file_name, NX, NE) : 
    e_cols = []
    for ei in range(NE) : 
        e_cols.append(str(ei))
    df = pd.read_csv(file_name, usecols=e_cols)
    var = df.values.tolist()
    return np.array(var)

def readInfo(file_number) : 
    if (file_number < 10) : 
        sfile_number = "0000"+str(file_number)
    if (file_number >= 10 and file_number < 100) : 
        sfile_number = "000"+str(file_number)
    if (file_number >= 100 and file_number < 1000) : 
        sfile_number = "00"+str(file_number)
    if (file_number >= 1000) : 
        sfile_number = "0"+str(file_number)
    
    file_name = "../data_out/"+"info"+"_"+sfile_number+".dat"
    df = open(file_name, "r").readlines()
    dic = {}
    for ii in range(len(df)) : 
        df[ii] = df[ii].replace("\n","").split(",")
        dic.update({df[ii][0]:df[ii][1]})
    return dic 

def readAxis(file_name) : 
    myfile = open(file_name,"r").readlines()
    for ii in range(len(myfile)) : 
        myfile[ii] = myfile[ii].replace("\n","")
        myfile[ii] = float(myfile[ii])
    return myfile


def getData(var, file_number) : 
    if (file_number < 10) : 
        sfile_number = "0000"+str(file_number)
    if (file_number >= 10 and file_number < 100) : 
        sfile_number = "000"+str(file_number)
    if (file_number >= 100 and file_number < 1000) : 
        sfile_number = "00"+str(file_number)
    if (file_number >= 1000) : 
        sfile_number = "0"+str(file_number)
    nx = fr.search("../parameters.dat", "NX")
    ne = fr.search("../parameters.dat", "NE")
    X = readAxis("../data_ini/X.dat")
    E = readAxis("../data_ini/E.dat")
    if (var == "Pcr" or var == "Ip" or var == "Im" or var == "Pe") : 
        data = readDataXE("../data_out/"+var+"_"+sfile_number+".dat", 2**int(nx), 2**int(ne))
    if (var == "Dcr") : 
        data = np.empty((len(X), len(E)))
        Ip = readDataXE("../data_out/"+"Ip"+"_"+sfile_number+".dat", 2**int(nx), 2**int(ne))
        Im = readDataXE("../data_out/"+"Im"+"_"+sfile_number+".dat", 2**int(nx), 2**int(ne))
        Db = readDataXE("../data_ini/"+"DBohm"+".dat", 2**int(nx), 2**int(ne))
        # for ii in range(len(Ip)) : 
        #     for jj in range(len(Ip[ii])) : 
        #         data[ii][jj] = Db[ii][jj]/(Ip[ii][jj] + Im[ii][jj])
        data = Db/(Ip + Im)
    return X, E, data



def show(variable, time, position, energy, 
         xlim = None, elim = None, 
         vlim = None, 
         source_center = 0, 
         fig_save = False, 
         info = False) : 
    
    X, E, data = getData(variable[0], 1)
    
    position_index = []
    energy_index   = []
    
    for ei in range(len(energy)) : 
        e_id = 0
        closest = np.inf 
        for ii in range(len(E)) : 
            if (abs(E[ii]-energy[ei]) < closest) : 
                e_id = ii
                closest = abs(E[ii]-energy[ei])
        energy_index.append(e_id)
    
    for xi in range(len(position)) : 
        x_id = 0
        closest = np.inf 
        for ii in range(len(X)) : 
            if (abs(X[ii]-position[xi]) < closest) : 
                x_id = ii
                closest = abs(X[ii]-position[xi])
        position_index.append(x_id)
    
    
    
    
    color = colorArray("blue", "red", n = len(time))
    
    
    
    
    size_x = 6
    size_y = 4
    sub_x  = 2
    sub_y  = max(len(position), len(energy))
    
    

    
    e_index = 1 
    x_index = 0 
    
    
    for vi in range(len(variable)) : 
        print ("We are preparing the variable : "+str(variable[vi]))
        fig = plt.figure(figsize=(size_x*sub_x,size_y*sub_y))
        gs = gridspec.GridSpec(ncols= sub_x, nrows = sub_y, figure = fig )
        gs.update(wspace=0.2, hspace=0.3) # set the spacing between axes. 
        for ei in range(len(energy)) : 
            print ("Spatial distribution subplot n°"+str(int(ei+1))+" over "+str(len(energy)))
            ax = fig.add_subplot(gs[ei, e_index])
            for ti in range(len(time)) : 
                X, E, data = getData(variable[vi], time[ti])
                if (source_center > 0) : 
                    ax.loglog(np.array(X)/cst.pc - source_center, data.T[energy_index[ei]], c = color[ti])
                    if (info) : 
                        loc_info = readInfo(time[ti])
                        # print ("Rsh = ",float(loc_info.get("R_sh"))/cst.pc," pc")
                        ax.axvline(float(loc_info.get("R_sh"))/cst.pc, c = color[ti])
                else : 
                    ax.semilogy(np.array(X)/cst.pc, data.T[energy_index[ei]], c = color[ti])
                    if (xlim) : ax.set_xlim(xlim[0], xlim[1])
                if (vlim) : ax.set_ylim(vlim[vi][0], vlim[vi][1])
                del data
            ax.set_title("Energy = "+str(round(E[energy_index[ei]]/cst.GeV,3))+" GeV")
             
        for xi in range(len(position)) : 
            print ("Energy spectra subplot n°"+str(int(xi+1))+" over "+str(len(position)))
            ax = fig.add_subplot(gs[xi, x_index])
            for ti in range(len(time)) : 
                X, E, data = getData(variable[vi], time[ti])
                ax.loglog(np.array(E)/cst.GeV, data[position_index[xi]], c = color[ti])
                if (elim) : ax.set_xlim(elim[0], elim[1])
                if (vlim) : ax.set_ylim(vlim[vi][0], vlim[vi][1])
                del data 
            ax.set_title("Position = "+str(round(X[position_index[xi]]/cst.pc,1))+" pc")
            
    
    

        
        if (fig_save) : 
            print ("Saving data as : "+variable[vi]+".pdf")
            plt.savefig(variable[vi]+".pdf")
            
        
        print ("========================================")








###############################################################################
# SHOW DATA                                                                   #
###############################################################################


# Experiments 
show(["Pcr", "Pe"], [0, 1, 2, 3, 4, 5, 6, 7], 
      [1040.*cst.pc, 1125.*cst.pc, 1300.*cst.pc, 1450.*cst.pc, 1600.*cst.pc, 1850.*cst.pc],
      [10.*cst.GeV, 100.*cst.GeV, 1.*cst.TeV, 3.*cst.TeV, 10.*cst.TeV, 20.*cst.TeV],
      elim = [1e0, 1e5],
      vlim = [[1e-26, 1e-8], [1e-26, 1e-8]],
      source_center = 1000.,
      fig_save = True)

show(["Dcr", "Ip", "Im"], [0, 1, 2, 3, 4, 5, 6, 7], 
      [1040.*cst.pc, 1125.*cst.pc, 1300.*cst.pc, 1450.*cst.pc, 1600.*cst.pc, 1850.*cst.pc],
      [10.*cst.GeV, 100.*cst.GeV, 1.*cst.TeV, 3.*cst.TeV, 10.*cst.TeV, 20.*cst.TeV],
      source_center = 1000.,
      fig_save = True)





# file_name = "../data_out/"+"info"+"_"+"00002"+".dat"
# df = open(file_name, "r").readlines()
# dic = {}
# for ii in range(len(df)) : 
#     df[ii] = df[ii].replace("\n","").split(",")
#     dic.update({df[ii][0]:df[ii][1]})
# NX = 2**12
# NE = 2**8

# data = readDataXE(file_name, NX, NE)

# data = readInfo(1)















