# -*- coding: utf-8 -*-
"""
Created on Wed May 31 05:34:10 2023

@author: SSI-User
"""
import numpy as np
import math
import matplotlib.pyplot as plt

d_N = np.array([4, 10, 30, 80, 280, 400, 1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000, 100000])
gradient_0 = np.array([-1, 3759, 4620, 4892, 5009, 5023, 5065, 5080, 5080,5080,5080,5080,5080,5080,5080])

fig = plt.figure(figsize=(10,5))

plt.plot(d_N, gradient_0,'o',markersize=10)

plt.xscale('log')
plt.minorticks_on()
plt.grid(b=True, which='major', color='black', linestyle='-')
plt.grid(b=True, which='minor', color='gray', linestyle='--')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title("Values of $\lambda$ for each iteration ",fontsize=25)
plt.xlabel('Number of discrete points (N)',fontsize=18)
plt.ylabel("Gradient at initial point ($\lambda$)",fontsize=18)
plt.savefig(r'C:\行星科学\Matlab\potential_data\Iteration_ide\number_gradient.png', dpi=300)
plt.show()