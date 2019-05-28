#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
FILE     : SETUP.PY 
FONCTION : Generate the folder of work, this folder will contain the initial values, 
           the properties of the simulation, the output data and the compiled source code  

"""

import os 
import sys 

# Input variables
DIMENSION = 2
NX        = 4
NY        = 2 
GRID      = "CARTESIAN"
OPENMP    = "ON"
NPROC     = 1 

# Other methods 




# Folders creation 
data_ini = "data_ini" 
data_out = "data_out"
params   = "params"
try : 
    os.mkdir(data_ini)
except : 
    pass
try : 
    os.mkdir(data_out)
except : 
    pass
try : 
    os.mkdir(params)
except : 
    pass 


# 