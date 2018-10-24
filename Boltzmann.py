# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:46:59 2018

@author: Ustun
"""
import math

def BoltzmannIntensity(Temperature,Energy):
    #returns normalized boltzmann relative population differences of two states
    boltzmannC = 8.6173303e-2 #boltzmann constant in meV/Kelvin
    I_1 = 1-1/(math.exp(-abs(Energy)/(boltzmannC*Temperature))+1)
    I_3 = 1/(math.exp(-abs(Energy)/(boltzmannC*Temperature))+1)
    return I_1,I_3