#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 09:51:47 2018

@author: lbrahimi
"""

import numpy as np
#
import matplotlib.pyplot as plt 
import pandas as pd
import os 

def search(file_name, variable) : 
    var = []
    file_ = open(file_name,"r").readlines()
    for ii in range(len(file_)) : 
        file_[ii] = file_[ii].split()
    for ii in range(len(file_)) : 
        if (variable == file_[ii][0]) : 
            var.append(file_[ii][1])
    return var[0]





def fileWrite(file_name, variables = {}, path='./', ext='.dat') : 
    myfile = open(path+file_name+ext,"w")
    
    for key, value in variables.items():
        myfile.write(str(key)+"\t"+str(value)+"\n")



def write1D(file_name, nx=None, ne=None, variable=None, path="./") : 
    data = {"Index": []}
    for ei in range(0, 2**ne, 1) : 
        temp_data = {ei : []}
        data.update(temp_data)
    df = pd.DataFrame(data)
    for xi in range(2**nx) : 
        loc_data = {}
        for ei in range(0, 2**ne, 1) : 
            temp_data = {"Index":xi, ei:variable[xi][ei]}
            loc_data.update(temp_data)
        df = df.append(loc_data, ignore_index=True)
    df.to_csv(path_or_buf=path+file_name)

def write2D(file_name, nx = None, ny = None, ne = None, variable = None, path = "./") : 
    import hilbert as h 
    h.hilbertData(nx, ny, ne, tensor = variable, file_name = file_name, folder_name = path)

def write1Daxis(file_name, variable=None, nx=None, path="./") : 
    myfile = open(path+file_name,"w")
    for ii in range(len(variable)) : 
        myfile.write(str(float(variable[ii]))+"\n")
    


# Example 
#nx = 3
#ne = 2
#
#tensor = np.zeros((2**nx, 2**ne)) 
#for xi in range(2**nx) : 
#    for ei in range(2**ne) : 
#        tensor[xi][ei] = ei#(xi+1)*(yi+1)*(ei+1)
#
#write1D("file.dat", nx=nx, ne=ne, variable=tensor, path="./")




