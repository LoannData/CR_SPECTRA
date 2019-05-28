#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:50:52 2019

@author: lbrahimi
"""
import pandas as pd
import numpy as np 


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
        "YIndex"       : [],
        "H(i+1)"       : [],
        "H(i-1)"       : [],
        "H(j+1)"       : [],
        "H(j-1)"       : []}
df = pd.DataFrame(data)


# We get our Hilbert curve 
hilbert = hi.getRHilbert2D(0., 0., NX, NY)

for ii in range(len(hilbert)) : 
    
    print "ii = ",ii,"/",len(hilbert)
    
    
    xi = int(hilbert[ii][0])
    yi = int(hilbert[ii][1])
    
    hup    = -1
    hdown  = -1
    hleft  = -1
    hright = -1
    for jj in range(len(hilbert)) : 

        
        xj = int(hilbert[jj][0])
        yj = int(hilbert[jj][1])
        # (i, j+1)
        if (xi == xj and yi + 1 == yj) : 
            hup    = jj
        # (i, j-1)
        if (xi == xj and yi - 1 == yj) : 
            hdown  = jj
        # (i+1, j)
        if (xi + 1 == xj and yi == yj) : 
            hright = jj
        # (i-1, j)
        if (xi - 1 == xj and yi == yj) : 
            hleft  = jj 
    
    df = df.append({"HilbertIndex": ii, 
                    "XIndex": xi, 
                    "YIndex": yi,
                    "H(i+1)": hright,
                    "H(i-1)": hleft,
                    "H(j+1)": hup,
                    "H(j-1)": hdown}, ignore_index=True) 


df.to_csv(path_or_buf="../../bin/"+file_name)


