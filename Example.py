#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:23:22 2020

@author: tobiaskohler
"""

import debye



Fe3O4_Fd3m = "/DebyeFunctionAnalysis/commit-2020-11-27/Fd-3m.cif" 
Fe3O4_P43212 = "/DebyeFunctionAnalysis/commit-2020-11-27/P4_32_12.cif" 
CeO2 =  "/DebyeFunctionAnalysis/commit-2020-11-27/CeO2_674b.cif"
CaF2 = "/DebyeFunctionAnalysis/commit-2020-11-27/CaF2.cif"       


dry_run = False
  
debye.Experiment(wavelength = 0.21148, 
                 occupancies = {'Fe1':0.98, 'Fe2': 0.39, 'Fe3':0.93, 'Fe4': 0.92}, 
                 gradient = None,
                 x_par = "Q",
                 x_min = 5.0, 
                 x_max = 10.0, 
                 x_step = 0.01,
                 partitions = 400, 
                 diameter = 12,
                 shape = 'Sphere',
                 cif_file = Fe3O4_P43212, 
                 filename = "Magh_10nm_APB_Q5to10", 
                 plot_crystal=False,
                 APB = True,
                 plot_results= "Q",
                 plot_2D_crystal = False,
                 output_path="/output/",
                 dry_run=dry_run)


debye.Experiment(wavelength = 0.21148, 
                  occupancies = {'Fe1':0.98, 'Fe2': 0.39, 'Fe3':0.93, 'Fe4': 0.92}, 
                  gradient = None,
                  x_par = "Q",
                  x_min = 5.0, 
                  x_max = 10.0, 
                  x_step = 0.01,
                  partitions = 400, 
                  diameter = 12,
                  shape = 'Sphere',
                  cif_file = Fe3O4_P43212, 
                  filename = "Magh_10nm_noAPB_Q5to10", 
                  plot_crystal=False,
                  APB = False,
                  plot_results= "Q",
                  plot_2D_crystal = False,
                  output_path="/output/",
                  dry_run=dry_run)

debye.Experiment(wavelength = 0.21148, 
                 occupancies = {'Fe1':0.98, 'Fe2': 0.39, 'Fe3':0.93, 'Fe4': 0.92}, 
                 gradient = None,
                 x_par = "Q",
                 x_min = 10.0, 
                 x_max = 15.0, 
                 x_step = 0.01,
                 partitions = 400, 
                 diameter = 12,
                 shape = 'Sphere',
                 cif_file = Fe3O4_P43212, 
                 filename = "Magh_10nm_APB_Q10to15", 
                 plot_crystal=False,
                 APB = True,
                 plot_results= "Q",
                 plot_2D_crystal = False,
                 output_path="/output/",
                 dry_run=dry_run)


debye.Experiment(wavelength = 0.21148, 
                  occupancies = {'Fe1':0.98, 'Fe2': 0.39, 'Fe3':0.93, 'Fe4': 0.92}, 
                  gradient = None,
                  x_par = "Q",
                  x_min = 10.0, 
                  x_max = 15.0, 
                  x_step = 0.01,
                  partitions = 400, 
                  diameter = 12,
                  shape = 'Sphere',
                  cif_file = Fe3O4_P43212, 
                  filename = "Magh_10nm_noAPB_Q10to15", 
                  plot_crystal=False,
                  APB = False,
                  plot_results= "Q",
                  plot_2D_crystal = False,
                  output_path="/output/",
                  dry_run=dry_run)

# debye.Experiment(wavelength = 0.21148, 
#                   occupancies = {'Fe(oct)':0.88, 'Fe(tet)': 0.97}, 
#                   gradient = None,
#                   x_par = "Q",
#                   x_min = 1.0, 
#                   x_max = 5.0, 
#                   x_step = 0.01,
#                   partitions = 30, 
#                   diameter = 6,
#                   shape = 'Sphere',
#                   cif_file = Fe3O4_Fd3m, 
#                   filename = "Magh_5nm_Fd3m_noAPB", 
#                   plot_crystal=False,
#                   APB = False,
#                   plot_results= "Q",
#                   plot_2D_crystal = False,
#                   output_path="/output/",
#                   dry_run=dry_run)

# debye.Experiment(wavelength = 0.21148, 
#                   occupancies = {'Fe(oct)':0.88, 'Fe(tet)': 0.97}, 
#                   gradient = None,
#                   x_par = "Q",
#                   x_min = 1.0, 
#                   x_max = 5.0, 
#                   x_step = 0.01,
#                   partitions = 30, 
#                   diameter = 6,
#                   shape = 'Sphere',
#                   cif_file = Fe3O4_Fd3m, 
#                   filename = "Magh_5nm_Fd3m_APB", 
#                   plot_crystal=False,
#                   APB = True,
#                   plot_results= "Q",
#                   plot_2D_crystal = False,
#                   output_path="/output/",
#                   dry_run=dry_run)

# debye.Experiment(wavelength = 0.21148, 
#                   occupancies = {}, 
#                   gradient = None,
#                   x_par = "Q",
#                   x_min = 1.0, 
#                   x_max = 5.0, 
#                   x_step = 0.01,
#                   partitions = 600, 
#                   diameter = 20,
#                   shape = 'Sphere',
#                   cif_file = CeO2, 
#                   filename = "CeO2_10nm", 
#                   plot_crystal=False,
#                   APB = False,
#                   plot_results= "Q",
#                   plot_2D_crystal = False,
#                   output_path="/output/",
#                   dry_run=dry_run)


        

            

