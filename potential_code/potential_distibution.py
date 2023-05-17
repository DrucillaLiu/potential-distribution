# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 00:47:30 2023

@author: SSI-User
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

jeb = 6.39 # s^-1m^-2
e = 1.6*10**-19 # C
e0 = 8.85*10**-12 # F/m
E = 250 # eV
m = 9.11*10**-31 #kg
E = 250 # eV
V0 = 0  # V
Vw = -153 # V 
a = 0.1  # neb0/new0
b = 25   # E/kTew
kTse = 3 # eV
kTew = E/b # eV
neb0 = 4.26*10**12 # m^-3  electron beam at sheath edge 
nse0 = 0.9*jeb/(e*(kTse/(2*3.14*m))**0.5)*np.exp(-e*Vw/kTse) # m^-3  secondary electrons at sheath edge
new0 = neb0/a      # m^-3  plasma electrons at sheath edge
ns = neb0+nse0+new0 # m^-3  ion number denstiy at sheath edge
Ws = ns/(4*(nse0/(2*kTse)+new0/(2*kTew)-neb0/(4*E))) # min ion energy at sheath edge 

h = 1e-6     # step of x
xmin = 0     # vertical distance from sample (m)
xmax = 50  # (m) 0.8 cm
x = np.arange(xmin, xmax+h, h) # create a distance scaler 
V = np.zeros(x.shape[0])       # create a potential scaler
dV = np.zeros(x.shape[0])      # electric field scaler
ign = np.zeros(x.shape[0])
V[0] = 0            # initial condition potential at x=0
dV0 = 1e-3          # electric field at x=0 (dV/dx at x=0) 
def integralN(v):        # density integral of electron beam 
    ig_n = (e*v)**2/e0*(1+nse0/(2*kTse)+new0/(2*kTew)-neb0/(4*E)-ns/(4*Ws))
    return ig_n
for i in np.arange(V.shape[0] - 1):
    dV[i] = (dV0**2+2*integralN(V[i]))**0.5  
    V[i+1] = V[i] - h*dV[i]
    ign[i] = integralN(V[i])
                       
plt.plot(x*100, V, 'orange')
plt.xlabel('distance (cm)')
plt.ylabel('potential (V)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.grid(True)   