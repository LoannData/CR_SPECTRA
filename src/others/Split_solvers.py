#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:48:23 2019

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 

# Thomas Algorithm (Trigonal matrix inversion O(N))
def TDMA(a,b,c,d):
    n = len(d)
    w= np.zeros(n-1,float)
    g= np.zeros(n, float)
    p = np.zeros(n,float)

    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]

    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p





# Generalized diffusion solver (implicit method)
def generalized_diffusion(u, X, D, dt, theta = 0.5) : 
    R = np.zeros((1, NX+1))
    
    a = np.zeros(NX+1)
    b = np.zeros(NX+1)
    c = np.zeros(NX+1)  
    
    for i in range(1, NX) : 
        dx = 0.5*(X[i+1] - X[i-1])
        F1 = theta*dt/dx**2
        alpha_m = 0.5*(D[i] + D[i-1])
        alpha_p = 0.5*(D[i] + D[i+1])
        
        an = -F1*alpha_m
        bn = F1*alpha_p + F1*alpha_m + 1 
        cn = -F1*alpha_p
        rn = u[i] + (1-theta)*dt/dx**2*(alpha_p*(u[i+1] - u[i]) - alpha_m*(u[i] - u[i-1]))
        
        a[i] = an
        b[i] = bn
        c[i] = cn
        R[0][i] = rn
        
    a = np.delete(a, 0)
    c = np.delete(c, -1)
    
    
    # Cas i = NX
    dx = X[1] - X[0]
    F1 = theta*dt/dx**2
    alpha_m = 0.5*(D[NX] + D[NX-1])
    alpha_p = alpha_m
    a[-1] = -F1*alpha_m
    b[-1] = F1*alpha_p + F1*alpha_m + 1
    
    # Cas i = 0 
    dx = X[1] - X[0]
    F1 = theta*dt/dx**2
    alpha_p = 0.5*(D[0] + D[1])
    alpha_m = alpha_p
    c[0]  = -F1*alpha_p
    b[0]  =  F1*alpha_p + F1*alpha_m + 1
    
    
    
    dx = 1.*(X[1] - X[0])
    u_new = TDMA(a,b,c,R[0])
    u_new[0] = 0.
    u_new[-1] = 0.
    
    return u_new 





NX = 1000
Xmin = -1
Xmax = 1
X = np.linspace(Xmin, Xmax, NX+1)


D = np.ones(len(X))*1e-2
D = D*(1+10*np.exp(-(X/0.2)**2))

u = np.ones(len(X))*0.
for ii in range(len(u)) : 
    if (abs(X[ii]) < 0.1) : 
        u[ii] = 1.


plt.plot(X, u, c="black")






t_ini = 0.
dt = 1e-1
t_max = 1.#2*dt
t = t_ini

u_new_0 = u
while (t <= t_max) : 
    u_old = u_new_0 
    u_new_0 = generalized_diffusion(u_old, X, D, dt, theta=0.5)
    t += dt 

t_ini = 0.
#dt = 1e-1
t_max = 1.#2*dt
t = t_ini

u_new_1 = u
while (t <= t_max) : 
    u_old = u_new_1 
    u_new_1 = generalized_diffusion(u_old, X, D, dt, theta=1.)
    t += dt 


plt.plot(X, u_new_0, c="blue")
plt.plot(X, u_new_1, c="red")




        
