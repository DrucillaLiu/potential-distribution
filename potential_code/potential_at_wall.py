# -*- coding: utf-8 -*-
"""
Created on Tue May 16 02:53:13 2023

@author: SSI-User
"""

import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import csv


# Fixed Parameters
e = 1.6*10**-19                     # C
m = 9.11*10**-31                    # kg
M = 1.67*10**-27                    # kg
kTse = 3                            # eV
kTpe = 10                           # eV
kTi = 10                            # eV
npe0 =  1e+18                       # m^-3  plasma electron density at sheath edge
vi0 = (2*e*kTi/M)**0.5              # m/s   ion velocity at sheath edge
pi = 3.14
zetape = ((kTpe*e)/(2*pi*m))**0.5   # j
zetase = ((kTse*e)/(2*pi*m))**0.5   # j
epsilonmax = 1
Emax = 300                          # eV

# Variables 
alpha = np.array([10**i for i in range(-2, 16, 1)])                 # neb0 / new0 [10^-2,10^-1,...,10^15]
beta = np.arange(20,85,5)                                           # Eeb0 / kTpe [10,15,20,25,30,...,80] 
neb = alpha * npe0                                                  # m^-3  electron beam density at sheath edge
Eeb = beta * kTpe                                                   # eV    electron beam energy at sheath edge

# Initialize arrays containing roots and errors
Uw0 = 0.5                                                           # alpha = 0.1, beta = 15, estimited starting root |Uw0| < Eeb0
potential = np.zeros([beta.shape[0], alpha.shape[0]])
error = np.zeros([beta.shape[0], alpha.shape[0]])

# Loop for finding roots by changing variables
for i in np.arange(beta.shape[0]):                                  # Eeb0 / kTpe [20,25,30,...,80]
    for j in np.arange(alpha.shape[0]):                             # neb0 / new0 [10^-2,10^-1,...,10^15]
        
        neb0 = alpha[j] * npe0
        Eeb0 = beta[i] * kTpe
        P = np.array([npe0*zetape*vi0, npe0*zetape, neb0*(2*e*Eeb0/m)**0.5*vi0, 2*neb0*zetase*(2*e*Eeb0/m)**0.5, 2*neb0*zetase*(2*e*Eeb0/m)**0.5-neb0*vi0-npe0*vi0])
        
        def func(x, P):
            A,B,C,D,H = P
            g1 = epsilonmax * ((Eeb0+x) / Emax) * np.exp(2-2*((Eeb0+x) / Emax)**0.5)
            g2 = epsilonmax * ((0.5*m*(zetape*np.exp(x/kTpe))**2) / Emax) * np.exp(2-2*((0.5*m*(zetape*np.exp(x/kTpe))**2) / Emax)**0.5)
            g3 = np.exp(-x/kTse)
            g4 = np.exp(x/kTpe)
            return A * g1 * g3 * g4 + B * (g1*g4-g4) + C * g2 * g3 + D * g1 - H
        
        root = fsolve(func, Uw0, args=P, full_output=1)
        
        Uw0 = root[0]
        fvec = root[1]['fvec']
        
        potential[i,j] = Uw0
        error[i,j] = fvec

# Save data to a csv file
density = ['10^16','10^17','10^18','10^19','10^20','10^21','10^22','10^23','10^24','10^25','10^26',
          '10^27','10^28','10^29','10^30','10^31','10^32','10^33']
#energy = ['200','250','300','350','400','450','500','550','600','650','700','750','800']
pd.DataFrame( potential ).to_csv( "potentials.csv", header = density, index = True )

data_dict = {
        '10^-2': potential[:,0], '10^-1': potential[:,1], '10^0': potential[:,2],'10^1': potential[:,3],
        '10^2': potential[:,4],'10^3': potential[:,5],'10^4': potential[:,6],'10^5': potential[:,7],
        '10^6': potential[:,8],'10^7': potential[:,9],'10^8': potential[:,10],'10^9': potential[:,11],
        '10^10': potential[:,12],'10^11': potential[:,13],'10^12': potential[:,14],
        '10^13': potential[:,15],'10^14': potential[:,16],'10^15': potential[:,17]
       } 