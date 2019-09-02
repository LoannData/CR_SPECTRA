#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 15:09:13 2018

@author: brahimi
"""
import numpy as np 

def search(file_name, variable) : 
    var = []
    file_ = open(file_name,"r").readlines()
    for ii in range(len(file_)) : 
        file_[ii] = file_[ii].split()
    for ii in range(len(file_)) : 
        if (variable == file_[ii][0]) : 
            var.append(file_[ii][1])
#    print var
    return var[0]

