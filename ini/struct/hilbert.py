#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:50:52 2019

@author: lbrahimi
"""
import pandas as pd


import hilbertRectangular2D as hi

# Simulation box size 
nx = 7
ny = 5

NX = 2**nx
NY = 2**ny

file_name = "grid2Dx"+str(nx)+"y"+str(ny)+".dat"

#HilbertIndex "\t" XIndex "\t" YIndex "\t" i-1 "\t" i+1 "\t" j-1 "\t" j+1"
data = {"HilbertIndex" : [],
        "XIndex"       : [], 
        "YIndex"       : []}
df = pd.DataFrame(data)


# We get our Hilbert curve 
hilbert = hi.getRHilbert2D(0., 0., NX, NY)

for ii in range(len(hilbert)) : 
    df = df.append({"HilbertIndex" : ii, "XIndex"       : int(hilbert[ii][0]), "YIndex"       : int(hilbert[ii][1])}, ignore_index=True) 


df.to_csv(path_or_buf="../../bin/"+file_name)


