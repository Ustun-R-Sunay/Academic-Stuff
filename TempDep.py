# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:58:35 2018

@author: Ustun
"""
import numpy as np
import Boltzmann
from matplotlib import pyplot as plt
import QuickHo
Dz = -5
Dx = 2.5E-4
LSc = -1
g_l = 0.1
H = QuickHo.CreateHamiltonianCF(Dz,Dx,LSc)
E_0,V_0 = QuickHo.GetStates(H)

gx1,gy1,gz1 = QuickHo.GetGValuesDeg(V_0[:,0],V_0[:,1],g_l)
gx2,gy2,gz2 = QuickHo.GetGValuesDeg(V_0[:,2],V_0[:,3],g_l)

#plt.plot(k,gx1,color = 'black',marker = 'o', linestyle = 'solid')
#plt.plot(k,gz1,color = 'blue',marker = 'o', linestyle = 'solid')

B_theta = np.linspace(0,180,1000)
gT1 = QuickHo.roadmap(gz1,gx1,B_theta*np.pi/180,0)
gT2 = QuickHo.roadmap(gz2,gx2,B_theta*np.pi/180,0)

plt.subplot(1,2,1)
plt.plot(B_theta,gT1,label = 'Band 1')
plt.plot(B_theta,gT2, label = 'Band 2')
plt.legend(loc = 9)
plt.xlabel("Angle (degrees)")
plt.ylabel("g-factor")
plt.title('$\Delta_z = {0}meV, \Delta_x = {1} meV,$ \n $\lambda = {2} meV, g_l = {3}$'.format(Dz,Dx,LSc,g_l))
#plt.axis([0,180,0,2.3])
Temp_range = np.linspace(0.1,10,50,endpoint='true')
I_1 = np.zeros(Temp_range.__len__())
I_3 = np.zeros(Temp_range.__len__())
c = 0
for i in Temp_range:
    I_1[c],I_3[c] = Boltzmann.BoltzmannIntensity(i,E_0[0]-E_0[2]) 
    c=c+1
x = 300
g_ave_d0 = I_1*gT1[x]+I_3*gT2[x]
print(B_theta[x])
plt.subplot(1,2,2)
gExp = g_ave_d0[17]-g_ave_d0
plt.plot(Temp_range,gT1[x]-g_ave_d0)
#plt.xlabel('Temperature (Kelvin)')
#plt.ylabel('g-factor')
##plt.show()
##print(E_0[A])
##print(V_0[:,A])
