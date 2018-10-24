# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 15:39:56 2018

@author: Ustun
"""

# -*- coding: utf-8 -*-

import sympy
import numpy as np

def SymDiag(H):
    A = H.eigenvects
    print(A)


def roadmap(gpar,gperp,angles,offset_angle):
    #plot angular dependence of EPR signal with axial symmetry
    #refer to any EPR book (e.g. weil,abragahm  for more details
    gTheta = np.sqrt(gpar*gpar*np.square(np.cos(angles+offset_angle)) +  
                     gperp*gperp*np.square(np.sin(angles+offset_angle)))
    return gTheta


def CreateHamiltonianCF(Dz,Dx,LSc):
    #create Hamiltonian 
    H = sympy.zeros(6,6)
    H[0,0] = -Dx/6 + Dz/3 + LSc/2
    H[2,0] = Dx/2
    H[1,1] = -Dx/6 + Dz/3 - LSc/2
    H[4,1] = LSc/np.sqrt(2)
    H[3,1] = Dx/2
    H[0,2] = Dx/2
    H[2,2] = -Dx/6 + Dz/3 - LSc/2
    H[5,2] = LSc/np.sqrt(2)
    H[1,3] = Dx/2
    H[3,3] = -Dx/6 + Dz/3 + LSc/2
    H[4,4] = Dx/3 -2*Dz/3
    H[1,4] = LSc/np.sqrt(2)
    H[2,5] = LSc/np.sqrt(2)
    H[5,5] = Dx/3 -2*Dz/3
    
    return H

def GetStates(H):
    #diagonalize Hamiltonian, Eigenvalues are "allowed energies" eigenvectors
    #describe the quantum state
    
    #some elements of H are very small and when diagonalization occurs,
    # you get incorrect eigenvectors. To fix this, multiply by a factor
    #E_0,V_0 = SymDiag(sympy.Matrix(H))
    factor = 1000
    H = H*factor
    
    E_0,V_0 = np.linalg.eig(H)
    #E_0,V_0 = A.eigenvects
    E_0 = E_0/factor
    #values aren't sorted by eigenvalues, do this manually
    A = np.argsort(E_0)
    E_0 = E_0[A]
    V_0 = V_0[:,A]
    
    return E_0,V_0

def GetGValuesDeg(wfa,wfb,gl):
    #calculates g-factors when ground state is orbitally degenerate.
    ge = 2.002319 #electron free g-factor value
    #The expressions can be simplified, but it is important to have physics be "readable"
    gx = 2*abs((wfb[0]*(gl*wfa[4]/np.sqrt(2)+ge*0.5*wfa[1]))+ wfb[1]*(gl*wfa[5]/np.sqrt(2)+ge*0.5*wfa[0]) + \
               wfb[2]*(gl*wfa[4]/np.sqrt(2)+ge*0.5*wfa[3]) + wfb[3]*(gl*wfa[5]/np.sqrt(2)+ge*0.5*wfa[2]) + \
               wfb[4]*(gl*(wfa[0]+wfa[2])/np.sqrt(2)+ge*0.5*wfa[5]) + wfb[5]*(gl*(wfa[1]+wfa[3])/np.sqrt(2)+ge*0.5*wfa[4]))
    
    gy = 2*abs(wfb[0]*(gl*wfa[4]/np.sqrt(2)+ge*0.5*wfa[1]))+ 2*abs(wfb[1]*(gl*wfa[5]/np.sqrt(2) - ge*0.5*wfa[0])) + \
               2*abs(wfb[2]*(-gl*wfa[4]/np.sqrt(2)+ge*0.5*wfa[3])) + 2*abs(wfb[3]*(-gl*wfa[5]/np.sqrt(2) + ge*0.5*wfa[2])) + \
               2*abs(wfb[4]*(gl*(wfa[0]-wfa[2])/np.sqrt(2)-ge*0.5*wfa[5])) + 2*abs(wfb[5]*(gl*(wfa[1]-wfa[3])/np.sqrt(2)+ge*0.5*wfa[4]))
    
    gz = 2*abs(wfb[0]*wfb[0]*(gl+ge/2)+wfb[1]*wfb[1]*(gl-ge/2)+wfb[2]*wfb[2]*(gl+ge/2)+wfb[3]*wfb[3]*(-gl-ge/2)+ \
               wfb[4]*wfb[4]*(-gl-ge/2)+wfb[5]*wfb[5]*(gl))
    return gx,gy,gz
    

