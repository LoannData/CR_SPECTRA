#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 14:38:57 2019

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 

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

data = readDataXE("../data_out/Pcr_0165.dat", 2**11, 2**5)
X = readAxis("../data_ini/X.dat")
E = readAxis("../data_ini/E.dat")


plt.figure(figsize=(6, 8))
#plt.loglog(E, data[500])
plt.loglog(np.array(X)/3.086e18, data.T[1])
plt.show()