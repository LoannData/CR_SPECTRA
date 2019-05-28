#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 09:51:47 2018

@author: lbrahimi
"""

import numpy as np
import matplotlib.pyplot as plt 

def search(file_name, variable) : 
    var = []
    file_ = open(file_name,"r").readlines()
    for ii in range(len(file_)) : 
        file_[ii] = file_[ii].split()
    for ii in range(len(file_)) : 
        if (variable == file_[ii][0]) : 
            var.append(file_[ii][1])
    return var[0]



    