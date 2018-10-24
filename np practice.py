# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 12:18:45 2018

@author: Ustun
"""

import sympy
import SympyEPR

A = sympy.Matrix([[3, -2,  4, -2], [5,  3, -3, -2], [5, -2,  2, -2], [5, -2, -3,  3]])
H = SympyEPR.CreateHamiltonianCF(-5,0,-1)
#B = list(A.eigenvals().keys())
Vals = H.eigenvects()

#C = [list(tup[2][0]) for tup in Vals ]
#
#D = [float(tup[0]) for tup in Vals]
print(Vals[0][2][1])
#
#print(Vals[0][2][0])  #returns eigenvalues and their algebraic multiplicity
#print(A.eigenvects())  #returns eigenvalues, eigenvects