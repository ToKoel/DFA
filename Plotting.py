#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 18:47:49 2020

@author: tobiaskohler
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def plot_3D(crystal):
    x,y,z, element_list = crystal.return_coordinate_list()

    colors = []
    for e in element_list:
        if e =="Fe":
            colors.append((0.2,0.2,0.0))
        else:
            colors.append((0.8,0.0,0.0))
    colors = np.array(colors)

    fig = plt.figure(figsize=plt.figaspect(3.0)*3.0)
    ax = fig.gca(projection='3d', proj_type='ortho')
 

    ax.scatter(x,y,z, color=colors)

    plt.show()