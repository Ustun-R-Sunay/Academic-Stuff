# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 12:18:45 2018

@author: Ustun
"""

import sympy
import SympyEPR
import numpy as np

H = SympyEPR.CreateHamiltonianCF(-5,0,-1)
#B = list(A.eigenvals().keys())
Vals = H.eigenvects()

#C = [list(tup[2][0]) for tup in Vals ]
#
#list comprehension
D = np.matrix([i[2][mult]/i[2][mult].norm() for i in Vals for mult in range(i[1])])
E = np.array([i[0] for i in Vals for mult in range (i[1])])
F = np.argsort(E)
E = E[F]
G = D[F]
#for i in Vals:
#    for multiplicity in range(i[1]):
#        E.append(float(i[0]))
#print(Vals[2])
print(G,'\n')
#D = [float(tup[0]) for tup in Vals]
#
#print(Vals[0][2][0])  #returns eigenvalues and their algebraic multiplicity
#print(A.eigenvects())  #returns eigenvalues, eigenvects