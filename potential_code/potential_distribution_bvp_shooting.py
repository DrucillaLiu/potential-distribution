# -*- coding: utf-8 -*-
"""
Created on Tue May 30 02:30:44 2023

@author: SSI-User
"""

import numpy as np
import math
import matplotlib.pyplot as plt

pi = 3.14
e = -1.6*10**-19             # C
e0 = 8.85*10**-12            # F/m
ne = 3.6*10**12              # m^-3  primary electron beam density
eta = 0.1                    # backscatterd electron emission yield
sigma = 0.9                  # secondary electron emission yield
E0 = 350                     # eV  primary electron beam energy
Us = -217                    # V   surface potential
kTse = 3                     # eV
A = e/e0*ne*(1+eta)
B = e/e0*ne*sigma*0.5*(4*pi*(E0/kTse))**0.5

## BVP
N = 100000
h = 0.2/N
x = np.linspace(0,0.2,N+1)
U1 = np.zeros(N+1)
U2 = np.zeros(N+1)
Z1 = np.zeros(N+1)
Z2 = np.zeros(N+1)

lambda_app = [5080]                  # lambda_0

U1[0] = 0
U2[0] = lambda_app[0]
Z1[0] = 0
Z2[0] = 1

tol = 0.0001
k = 0
fig = plt.figure(figsize=(10,5))
while  k<10:
    k = k+1
    for i in range (0,N):
        U1[i+1] = U1[i]+h*U2[i]
        U2[i+1] = U2[i]+h*(A*(E0/(E0+U1[i]))**0.5+B*np.exp((Us-U1[i])/kTse))
    
        Z1[i+1] = Z1[i]+h*(Z2[i])
        Z2[i+1] = Z2[i]+h*((-0.5*A*E0**0.5*(E0+U1[i])**-1.5+B/kTse*np.exp((Us-U1[i])/kTse))*Z1[i])

    lambda_app.append(lambda_app[k-1]-(U1[N]-Us)/Z1[N])
    #plt.plot(x*100,U1,':o',label=r"$\lambda$ %s"%(k),markersize=5,linewidth=1)
    U2[0] = lambda_app[k]
    if abs(lambda_app[k]-lambda_app[k-1])<tol:
        break

#plt.plot(x*100,U1,'b:o', linewidth=2, markersize=5)  
plt.plot(lambda_app,'o',markersize=10)
plt.grid(True)
#plt.legend(loc='best')
#plt.minorticks_on()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.yticks([-200,-150,-100,-50,0, 50,100,150,200], fontsize=15)
plt.xlabel('Iterations k',fontsize=18)
plt.ylabel("Gradient at initial point $\lambda$",fontsize=18)
plt.title("Values of $\lambda$ for each iteration ",fontsize=25)
#plt.xlabel('Distance from initial point (cm)',fontsize=18)
#plt.ylabel('Potential (V)',fontsize=18)
#plt.title("Numerical Approximation of the potential distribution",fontsize=25) 
#plt.title(r"Iterations of $\lambda_k$",fontsize=25)
#plt.savefig('C:\行星科学\Matlab\potential_data\lambda_ide\lambda_-1_4.png', dpi=300) # _lambda0_N
#plt.savefig(r'C:\行星科学\Matlab\potential_data\numerticalpotential_ide\numpot_5080_100000.png', dpi=300)
plt.savefig(r'C:\行星科学\Matlab\potential_data\Iteration_ide\Iteration_5080_100000.png', dpi=300)
plt.show()









