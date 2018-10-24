# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 13:06:09 2018

@author: Ustun
"""

#application of QuickHo and Boltzmann code
import scipy as sp
import scipy.linalg

import numpy as np
a = np.array([[-800.21,-600.00],[-600.00,-1000.48]])
eigvals, eigvecs = sp.linalg.eig(a)

#define variables
Dz = -5  #axial crystal field in meV
Dx = 0   #Basal crystal field in meV
LSc = -1  #spin-orbit coupling strength in meV
g_l = 0.1 #covalency or "orbital reduction" factor, no units
