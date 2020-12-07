#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 18:47:49 2020

@author: tobiaskohler
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_3D(crystal=None, coordinates = None, element_list = None):
    if crystal != None:
        coordinates, element_list = crystal.return_coordinate_list()
    else:
        coordinates = coordinates
        element_list = element_list

    colors = []
    for e in element_list:
        if e =="Fe":
            colors.append((0.2,0.2,0.0))
        elif e == "O":
            colors.append((0.8,0.0,0.0))
        else:
            colors.append((0.0,0.0,0.0))
    colors = np.array(colors)

    fig = plt.figure(figsize=plt.figaspect(3.0)*3.0)
    ax = fig.gca(projection='3d', proj_type='ortho')
 
    x,y,z = zip(*coordinates)
    ax.scatter(x,y,z, color=colors)

    plt.show()

def plot_2D(crystal=None, coordinates = None, element_list = None):
    if crystal != None:
        coordinates, element_list = crystal.return_coordinate_list()
    else:
        coordinates = coordinates
        element_list = element_list

    colors = []
    for e in element_list:
        if e =="Fe":
            colors.append((0.2,0.2,0.0,0.5))
        elif e == "O":
            colors.append((0.8,0.0,0.0,0.0))
        else:
            colors.append((0.0,0.0,0.0))
    colors = np.array(colors)

    fig,ax = plt.subplots(1,1)
 
    x,y,z = zip(*coordinates)
    indices = []
    lim = crystal.diameter/2
    for n,i in enumerate(z):
        if i  >= (crystal.diameter/2)-lim and i <= (crystal.diameter/2)+lim:
            indices.append(n)
    xs = []
    ys = []
    cs = []
    
    for i in indices:
        xs.append(x[i])
        ys.append(y[i])
        cs.append(colors[i])
        
    ax.scatter(xs,ys, color=cs)
    ax.set_aspect('equal')
    plt.show()    


def plot(data,xlabel,ylabel):
    fig,ax = plt.subplots(1,1)
    for d in data.keys():
        x,y,linewidth, marker = data[d]
        ax.plot(x,y,linewidth=linewidth, marker=marker,label=d)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    
    
def plot_results(x_par,x,y):
    fig,ax = plt.subplots(1,1)
    ax.plot(x,y)
    if x_par == "Q":
        ax.set_xlabel("$Q [\AA^{-1}]$")
    elif x_par == "2Theta":
        ax.set_xlabel("$2\Theta$ [°]")
    ax.set_ylabel("Intensity [a.u.]")
    plt.show()