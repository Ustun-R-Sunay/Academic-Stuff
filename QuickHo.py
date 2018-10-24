# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:04:28 2018

@author: Ustun
"""
import sympy
import numpy as np
from matplotlib import pyplot as plt

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
    H = np.zeros((6,6))
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
    

#Example code
#    #clean this code up for publication...
#Dz = -1800
#LSc = -15
#g_l = 0.1
#gx1 = np.zeros(50)
#gy1 = np.zeros(50)
#gz1 = np.zeros(50)
#
#gx2 = np.zeros(50)
#gy2 = np.zeros(50)
#gz2 = np.zeros(50)
#c = 0 
#k = np.logspace(-3,2,50,endpoint=True)
##k = range(51)
#for i in k:
#    H = CreateHamiltonianCF(Dz,i,LSc)
#    E_0,V_0 = GetStates(H)
#    gx1[c],gy1[c],gz1[c] = GetGValuesDeg(V_0[:,0],V_0[:,1],0.1)
#    gx2[c],gy2[c],gz2[c] = GetGValuesDeg(V_0[:,2],V_0[:,3],0.1)
#
#    c = c+1
#plt.plot(k,gx1,color = 'black',marker = 'o', linestyle = 'solid')
#plt.plot(k,gz1,color = 'blue',marker = 'o', linestyle = 'solid')
#B_theta = np.linspace(0,np.pi,1000)
#gT1 = roadmap(gz1[0],gx1[0],B_theta,0)
#gT2 = roadmap(gz2[0],gx2[0],B_theta,0)
#
#plt.subplot(1,2,1)
#plt.plot(B_theta*180/np.pi,gT1,label = 'Band 1')
#plt.plot(B_theta*180/np.pi,gT2, label = 'Band 2')
#plt.legend(loc = 9)
#plt.xlabel("Angle (degrees)")
#plt.ylabel("g-factor")
#plt.title('$\Delta_z = {0}meV, \Delta_x = {1} meV,$ \n $\lambda = {2} meV, g_l = {3}$'.format(Dz,k[0],LSc,g_l))
#plt.axis([0,180,0,2.2])
#Temp_range = np.linspace(0.1,10,50,endpoint='true')
#I_1 = np.zeros(Temp_range.__len__())
#I_3 = np.zeros(Temp_range.__len__())
#c = 0
#for i in Temp_range:
#    I_1[c],I_3[c] = Boltz.BoltzmannIntensity(i,E_0[0]-E_0[2]) 
#    c=c+1
#g_ave_d0 = I_1*gT1[0]+I_3*gT2[0]
#
#
#plt.subplot(1,2,2)
#plt.plot(Temp_range,g_ave_d0)
#plt.xlabel('Temperature (Kelvin)')
#plt.ylabel('g-factor')
##plt.show()
##print(E_0[A])
##print(V_0[:,A])