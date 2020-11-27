#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:23:22 2020

@author: tobiaskohler
"""

import XtlParser
import Crystal
import DFA
import numpy as np
from Plotting import plot_3D

import matplotlib.pyplot as plt

# atom = Crystal.Atom(coordinates = (0.1,0.1,0.1))
# atom2 = Crystal.Atom(coordinates = (0.3,0.2,0.4))

# atom.calculate_distance(atom2)
# print(atom.distances)




filename ="Fd-3m.xtl"
unitcell = XtlParser.XtlParser(filename)

crystal = Crystal.Crystal(diameter=5, unitcell = unitcell.coordinatelist, 
                          elements= unitcell.elementlist, latticepar = unitcell.latticepar)

# print(len(crystal.atoms))

crystal.generate_APB()
crystal.cut_sphere()
crystal.calculate_formfactor_matrix()
crystal.calculate_distance_matrix()


calc = DFA.DSE_calculator(crystal,3.5,14,0.2,0.2115,8.49)
q,theta,I = calc.calculateDSE()

plt.plot(theta,I)
plt.show()

# print(len(crystal.atoms))

#plot_3D(crystal)