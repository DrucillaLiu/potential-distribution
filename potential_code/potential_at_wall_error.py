# -*- coding: utf-8 -*-
"""
Created on Wed May 17 20:47:47 2023

@author: SSI-User
"""
# Import packages
from dash import Dash, html, dash_table
import pandas as pd
import numpy as np

# Incorporate data
df = pd.read_csv('https://raw.githubusercontent.com/DrucillaLiu/potential-distribution/main/potential_data/potentials_at_wall.csv')

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

# Loop for recording errors
potential = np.zeros([beta.shape[0], alpha.shape[0]])
error = np.zeros([beta.shape[0], alpha.shape[0]])
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
        
        potential[i,j] = df.iloc[i,j+1]
        error[i,j] = func(potential[i,j], P)
        
# Save data to a csv file
density = ['neb=10^16','neb=10^17','neb=10^18','neb=10^19','neb=10^20','neb=10^21','neb=10^22','neb=10^23','neb=10^24','neb=10^25','neb=10^26',
          'neb=10^27','neb=10^28','neb=10^29','neb=10^30','neb=10^31','neb=10^32','neb=10^33']
energy = ['Eeb=200','Eeb=250','Eeb=300','Eeb=350','Eeb=400','Eeb=450','Eeb=500','Eeb=550','Eeb=600','Eeb=650','Eeb=700','Eeb=750','Eeb=800']
df_error = pd.DataFrame(error, index=energy, columns=density)
df_p = pd.DataFrame(potential, index=energy, columns=density)
pd.DataFrame(error).to_csv( "potential_at_wall_error.csv", header = density, index = True )
pd.DataFrame(df_p).to_csv( "potential_at_wall.csv", index = True)
