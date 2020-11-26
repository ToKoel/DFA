#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:23:22 2020

@author: tobiaskohler
"""

import XtlParser
import Crystal
import numpy as np
from Plotting import plot_3D




filename ="/Users/tobiaskohler/PhD/Crystal structures/P4_332_unit_cell.xtl"
unitcell = XtlParser.XtlParser(filename)

crystal = Crystal.Crystal(diameter=4, unitcell = unitcell.coordinatelist, 
                          elements= unitcell.elementlist)

print(len(crystal.atoms))

crystal.generate_APB()
crystal.cut_sphere()

print(len(crystal.atoms))

plot_3D(crystal)