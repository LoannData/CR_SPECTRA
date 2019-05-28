#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:46:13 2019

@author: lbrahimi
"""

import pandas as pd
import numpy as np 

file_name = "../../bin/grid2Dx3y2.dat"
data = pd.read_csv(file_name, usecols=["HilbertIndex", "XIndex", "YIndex", "H(i+1)", "H(i-1)", "H(j+1)", "H(j-1)"])



#for ii in range(len(data)) : 
#    xi = data["XIndex"].iloc[ii]
#    yi = data["YIndex"].iloc[ii]
#    hi = data["HilbertIndex"].iloc[ii]
#    
#    xip = data["XIndex"].
    
    
    