#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 19:02:27 2020

@author: tobiaskohler
"""
import numpy as np

#Cromer-Mann coefficients for atomic formfactors
Fe3_coeff = [10.0333,4.36001,7.90625,0.26250,
             4.20562,9.35847,0.55048,20.4105,
             0.30429]
Fe2_coeff = [10.1270,4.44133,7.78007,0.27418,
             4.71825,10.1451,0.89547,24.8302,
             0.47888]
Fe_coeff = [11.769,7.357, 3.522, 2.305,
              4.761,0.307,15.354, 76.881,
              1.037]
O2_coeff = [3.7504,16.5131,2.9429,6.5920,
            1.5430,0.3192,1.6209,43.3486,
            0.2421]

def formfactor(theta_min, theta_max,step, wavelength, element):
    theta = np.arange(theta_min, theta_max, step)*np.pi/180
    s = np.sin(theta)/wavelength
    coeff = None
    if element == "Fe3+":
        coeff = Fe3_coeff
    elif element == "Fe2+":
        coeff = Fe2_coeff
    elif element =="O":
        coeff = O2_coeff
        
    f = (coeff[0]*np.exp(-coeff[4]*s**2)+
        coeff[1]*np.exp(-coeff[5]*s**2)+
        coeff[2]*np.exp(-coeff[6]*s**2)+
        coeff[3]*np.exp(-coeff[7]*s**2)+Fe3_coeff[8])
    return np.array(f)


 

class DSE_calculator():
    def __init__(self, crystal, t_min, t_max,step, wavelength, latticepar):
        self.crystal = crystal
        self.t_min = t_min
        self.t_max = t_max
        self.step = step
        self.wavelength = wavelength
        self.latticepar = latticepar
        self.n_fe = 0
        self.n_o = 0
        
        formfactor_Fe3 = formfactor(t_min, t_max,step, wavelength, "Fe3+")
        formfactor_Fe2 = formfactor(t_min, t_max,step, wavelength, "Fe2+")
        formfactor_O = formfactor(t_min, t_max,step, wavelength, "O")
        
        self.formfactor_Fe3_sq = formfactor_Fe3**2
        self.formfactor_O_sq = formfactor_O**2
        self.formfactor_Fe3_O = formfactor_Fe3*formfactor_O
        
        # Extract distance values
        # extract diagonal values, all are zero
        self.rii = np.diagonal(self.crystal.distance_matrix)
        #extract lower half of matrix without diagonal (k=-1)
        self.rij = np.tril(self.crystal.distance_matrix, k=-1)
        # flatten matrix into array
        self.rij = self.rij.flatten()
        # remove zero entries
        self.rij = self.rij[self.rij != 0]
        # multiply by lattice parameter to obtain real distances
        self.rij *= self.latticepar
        
        # Extract formfactor values
        # extract diagonal values of matrix
        self.fii = np.diagonal(self.crystal.formfactor_matrix) 
        
        # extract lower half of matrix without diagonal (k=-1)
        self.fij = np.tril(self.crystal.formfactor_matrix, k=-1)
        # flatten matrix into array
        self.fij = self.fij.flatten()
        # remove zero entries
        self.fij = self.fij[self.fij != 0]
        # replace placeholders with formfactors
        ffs = []
        for i in self.fij:
            if i == 4:
                ffs.append(self.formfactor_Fe3_sq)
            if i == 6:
                ffs.append(self.formfactor_Fe3_O)
            if i == 9:
                ffs.append(self.formfactor_O_sq)
        self.fij = np.array(ffs)
        
        # replace placeholders with formfactors
        for i in self.fii:
            if i == 4:
                self.n_fe += 1
            if i == 9:
                self.n_o += 1


        
    def calculateDSE(self):
        
        # generate q-array
        theta = np.arange(self.t_min, self.t_max, self.step)
        q = np.sin(theta*np.pi/180)/self.wavelength
        
        # convert distances array into vertical array and multiply with q
        rij_reshaped = self.rij.reshape(len(self.rij),1)
        # qr contains: [[qr11], [qr12], ... , [qrij]]
        qr = q*rij_reshaped
        
        I_diff = 0
        for a,n in enumerate(qr):
            I_diff += self.fij[a]*np.sin(n)/n
            
        I_same = (self.n_fe*self.formfactor_Fe3_sq
                  + self.n_o*self.formfactor_O_sq)
  
        
        I = 2*I_diff + I_same

        return q,theta,I

        
