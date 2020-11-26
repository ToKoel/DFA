#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:32:09 2020

@author: tobiaskohler
"""
import numpy as np
from itertools import product

class Atom():
    
    def __init__(self, coordinates=(0.0,0.0,0.0), element="Fe"):
        self.coordinates = np.array(coordinates)
        self.element = element
        
        
class Crystal():
    def __init__(self, diameter, unitcell, elements):
        self.diameter = diameter
        self.radius = diameter/2.0
        self.atoms = list()
        self.unitcell = unitcell
        self.elements = elements
        
        for a,e in zip(self.unitcell, self.elements):
            for f in product(range(self.diameter), range(self.diameter), range(self.diameter)):
                self.atoms.append(Atom(coordinates= a + np.array(f), 
                                            element = e))
        self.atoms = np.array(self.atoms)
        
    def generate_APB(self):
        for atom in self.atoms:
            atom.coordinates -= self.radius
            if (atom.coordinates[0]-atom.coordinates[1]) > 0.0:
                atom.coordinates += np.array([0.25, 0.25, 0.0])
            atom.coordinates += self.radius
            
    def cut_sphere(self):
        cut_crystal = []
        for atom in self.atoms:
            atom.coordinates -= self.radius
            if (np.dot(atom.coordinates, atom.coordinates) < self.radius**2):
                atom.coordinates += self.radius
                cut_crystal.append(atom)
        self.atoms = cut_crystal
        
    def return_coordinate_list(self):
        x = []
        y = []
        z = []
        element_list = []
        for atom in self.atoms:
            x.append(atom.coordinates[0])
            y.append(atom.coordinates[1])
            z.append(atom.coordinates[2])
            element_list.append(atom.element)         
        return x,y,z, element_list
            
            
        
                
    def rotate(self,alpha,beta,gamma):
        alpha = np.pi/180.0
        beta = np.pi/180.0
        gamma = np.pi/180.0
        
        rotation_x = np.matrix([1.0,0.0,0.0],
                               [0.0, np.cos(alpha), -np.sin(alpha)],
                               [0.0, np.sin(alpha), np.cos(alpha)])
        
        rotation_y = np.matrix([np.cos(beta), 0.0, np.sin(beta)],
                               [0.0, 1.0, 0.0],
                               [-np.sin(beta), 0.0, np.cos(beta)])
        
        rotation_z = np.matrix([np.cos(gamma),-np.sin(gamma),0.0],
                               [np.sin(gamma), np.cos(gamma), 0.0],
                               [0.0, 0.0,1.0])
        
        for atom in self.atoms:
            atom.coordinates -= self.radius
            atom.coordinates = np.matmul(rotation_x, atom.coordinates)
            atom.coordinates = np.matmul(rotation_y, atom.coordinates)
            atom.coordinates = np.matmul(rotation_z, atom.coordinates)
            atom.coordinates += self.radius
        