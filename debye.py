#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:59:36 2020

@author: tobiaskohler
"""

import crystal
import dse
import plotting

class Experiment():
    def __init__(self,wavelength,
                
                 x_step,
                 x_min,
                 x_max,
                 partitions,
                 diameter,
                 unitcell,
                 filename,
                 
                 x_par = "Theta",
                 occupancies = {},
                 shape="Sphere",
                 APB=False,
                 gradient=None,
                 plot_crystal=False,
                 plot_results=None,
                 output_path = "output/",
                 plot_2D_crystal = False):
  
        if occupancies == {}:
            occ = False
        else:
            occ = True
            
        print("===============================")
        print("Experiment started")
        print("===============================")
        
        Crystal = crystal.Crystal(diameter, unitcell, occupancies, shape)
        
        Crystal.build_nanoparticle(APB = APB,
                                   occ = occ,
                                   gradient = gradient,
                                   plot = plot_crystal)
    
        if plot_2D_crystal:
            plotting.plot_2D(Crystal)
        
        DSE_calculator = dse.DseCalculator(Crystal,x_par, x_min, x_max, 
                                            x_step, wavelength, 
                                            partitions, output_path)
        DSE_calculator.perform_calculation(filename = filename, plot = plot_results)
        
        print("===============================")
        print("Experiment finished")
        print("===============================")