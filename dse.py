#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 19:02:27 2020

@author: tobiaskohler

Debye scattering equation (Dse) module

Generate atomic formfactors and calculate the Debye scattering equation for 
nanocrystals.

Classes:
    
    AtomicFormfactorCalculator
    DseCalculator
    
Variables:
    Fe3p_coeff
    Fe2p_coeff
    Fe_coeff
    Ce4p_coeff
    O2m_coeff
    F1m_coeff
    Ca2p_coeff 
"""
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import numexpr as ne
import plotting

#Cromer-Mann coefficients for atomic formfactor from 
# https://www.classe.cornell.edu/~dms79/x-rays/f0_CromerMann.txt

Fe3p_coeff = [11.17640,  7.386300 ,  3.394800 ,  7.2400004E-02 , 
             4.614700 ,  0.3005000 ,  11.67290 , 38.55660,
              0.9707000]
Fe2p_coeff = [11.04240,  7.374000 ,  4.134600 ,  0.4399000,
             4.653800 ,  0.3053000  , 12.05460 , 31.28090,
              1.009700]
Fe_coeff = [11.76950 ,  7.357300 ,  3.522200 ,  2.304500 ,
             4.761100  , 0.3072000 ,  15.35350 , 76.88050,
             1.036900]

Ce4p_coeff = [20.32350, 19.81860, 12.12330, 0.1445830, 
             2.659410, 0.2188500, 15.79920, 62.23550,
             1.591800]

O2m_coeff = [4.758000,  3.637000,   0.000000,   0.000000, 
             7.831000,   30.05000,   0.000000,  0.000000,
             1.594000]

F1m_coeff = [3.632200 ,  3.510570,   1.260640 ,  0.9407060 , 
             5.277560 ,  14.73530 ,  0.4422580 , 47.34370,
              0.6533960]

Ca2p_coeff = [15.63480  , 7.951800 ,  8.437200  , 0.8537000 ,
                -7.4000000E-03 ,  0.6089000  , 10.31160 , 25.99050,
                -14.87500]

class AtomicFormfactorCalculator():
    """
    Class to generate formfactors.
    
    Attributes
    ----------
    formfactors : defaultdict(list)
        formfactors stored in a dictionary, can be retrieved with the element key.
    theta : ndarray
        Theta angles in degrees between min and max values and step size specified 
        in input.
    theta_rad : ndarray
        Theta angles in rad.
    q : ndarray
        Scattering vector q defined as 4pi*sin(theta)/lambda.
    q_4pi : ndarray
        Scattering vector q defined as sin(theta)/lambda.
    elements : list(str)
        List of strings representing the elements/ions to be calculated.
        
    Methods
    -------
    formfactor(q,element):
        Calculate the formfactor in range q of element.
    get_formfactor(element):
        Return formfactor array of element.
    get_theta():
        Return theta array.
    get_q():
        Return q array.
    plot_formfactors():
        Plot all formfactors in the range they were calculated.

    """ 
    def __init__(self,x_par, x_min, x_max, x_step, wavelength, elements):
        """
        Initialize the AtomicFormfactorCalculator with a specified q or Theta
        range.

        Parameters
        ----------
        x_par : str
            String label to set which x-values to use (2Theta or Q).
        x_min : float
            Lower limit for x.
        x_max : float
            Upper limit for x.
        x_step : float
            Step size of x..
        wavelength : float
            X-ray wavelength in Angstrom.
        elements : list of str
            List containing str labels of the elements to be calculated.

        Returns
        -------
        None.

        """
        self.formfactors = defaultdict(list)
        if x_par == "Theta":
            self.theta = np.arange(x_min, x_max, x_step)
            self.theta_rad = np.arange(x_min, x_max, x_step)*np.pi/180
            self.q = 4*np.pi*np.sin(self.theta*np.pi/180)/wavelength 
            self.q_4pi = np.sin(self.theta*np.pi/180)/wavelength
        elif x_par == "Q":
            self.q = np.arange(x_min, x_max, x_step)
            self.q_4pi = self.q/(4*np.pi)
            self.theta_rad = np.arcsin(self.q*wavelength/(4*np.pi))
            self.theta = self.theta_rad * 180/np.pi
        self.elements = elements
        for i in self.elements:
            self.formfactors[i].append(self.formfactor(self.q_4pi, i))

            
    def formfactor(self, q, element):
        """
        Calculate the formfactor in range q of element.
        
        The non-dispersive part of the atomic scattering factor is calculated 
        in the specified q range according to 
        f(k) = c+[SUM a_i*EXP(-b_i*(k^2))]; i=1,4
        with k = sin(theta)/lambda and c,a_i and b_i the Cromer-Mann coefficients
        given in the International Tables of Crystallography vol. 4 or vol C 
        (pg 500-502) and tabulated in 
        https://www.classe.cornell.edu/~dms79/x-rays/f0_CromerMann.txt. This 
        parametrization is only good up to sin(theta)/lambda < 2.0 [Angstrom^-1].
        For lambda = 0.21148 this means an upper limit of theta of 25 [°]. The
        correct parameters are chosen according to the element/ion name provided
        in the function argument.
        

        Parameters
        ----------
        q : ndarray
            Q array used for the calculation.
        element : str
            Element/Ion to be calculated.

        Returns
        -------
        ndarray
            Formfactor of element in the specified q range.

        """
        coeff = None
        if element == "Fe3+":
            coeff = Fe3p_coeff
        
        elif element == "Fe2+":
            coeff = Fe2p_coeff
        
        elif element =="O2-":
            coeff = O2m_coeff
            
        elif element =="Ce4+":
            coeff = Ce4p_coeff
            
        elif element =="F1-":
            coeff = F1m_coeff
        
        elif element =="Ca2+":
            coeff = Ca2p_coeff

    
        f = (coeff[0]*np.exp(-coeff[4]*q**2)+
        coeff[1]*np.exp(-coeff[5]*q**2)+
        coeff[2]*np.exp(-coeff[6]*q**2)+
        coeff[3]*np.exp(-coeff[7]*q**2)+coeff[8])
        return np.array(f)
    
    def get_formfactor(self, element):
        return np.array(self.formfactors[element][0])
    
    def get_theta(self):
        return self.theta
    
    def get_q(self):
        return self.q
    
    def plot_formfactors(self):
        fig,ax = plt.subplots(1,1)
        
        for i in self.elements:
            ax.plot(self.theta*2, self.get_formfactor(i), label=i)
  
        ax.legend()
        ax.set_xlim(left=0)
        ax.set_xlabel("$2\Theta (°)$")
        ax.set_ylabel("f (electrons)")
        plt.show()
            
class DseCalculator():
    """
    Class to calculate the Debye scattering equation.
    
    Attributes
    ----------
    distances : ndarray
        Array containing the pair wise distances.    
    f_ij : chararray
        Array containing the labels corresponding to atomic scattering factor
        products in the same order as the distances. (E.g. "FeO")     
    f_ii : chararray
        Array containing the diagonal elements of the atomic scattering factor
        matrix, i.e. the self correlations.      
    latticepar : float
        Lattice parameter a of the unit cell.      
    atom_numbers : dictionary
        Dictionary containing the total number of atoms grouped by elements.    
    t_min : float
        Lower limit of the theta range.      
    t_max : float
        Upper limit of the theta range.      
    step : float
        Theta step.      
    wavelength : float
        X-ray wavelength in Angstrom.      
    num_partitions : int
        Number of partitions to be performed prior to calculation.     
    intensity : ndarray
        Array containing the I(theta).
       
    Methods
    -------
    partition_crystal():
        Divide the crystal into parts to keep memory usage within RAM limits.
    calculate_Dse(distances, f_ij):
        Calculate the Debye scattering equation for non-equal pairs of atoms.
    calculate_Dse_total():
        Calculate the scattered intensity from all parts, including the self
        correlations.
    perform_calculation(filename=None, plot=True):
        Wrapper function to ensure results are saved and plotted or intentionally 
        discarded. 
    save_to_text(filename):
        Utility function to save q,theta and I to text file. 
        
    """ 
    
    def __init__(self,crystal,x_par,x_min, x_max,x_step, wavelength,num_partitions,output_path):
        """
        Initialize the DseCalculator.

        Parameters
        ----------
        crystal : Crystal object
            Crystal object obtained from crystal module.
        x_par : str
            Parameter to specify the x-values to be used (Theta or Q).
        x_min : float
            Lower limit of x range.
        x_max : TYPE
            Upper limit of x range.
        x_step : float
            Step size of x.
        wavelength : float
            X-ray wavelength in Angstrom.
        num_partitions : int
            Number of partitions to be performed prior to calculation.
        output_path: str
            Path for output files, only existing directories are accepted.

        Returns
        -------
        None.

        """
        self.cif_filename = crystal.cif_filename
        self.occupancies = crystal.occupancies
        self.diameter = crystal.diameter
        self.shape = crystal.shape
        self.distances, self.f_ij, self.f_ii = crystal.get_properties()
        self.latticepar = crystal.lattice_a 
        self.atom_numbers = crystal.atom_numbers
        self.x_par = x_par
        self.x_min = x_min
        self.x_max = x_max
        self.x_step = x_step
        self.wavelength = wavelength
        self.num_partitions = num_partitions
        self.intensity = 0
        self.output_path = output_path
       
        elements = ["O2-", "Fe3+", "Fe2+", "Ce4+", "Ca2+", "F1-"]

        formfactors = AtomicFormfactorCalculator(self.x_par,self.x_min, 
                                                   self.x_max,
                                                   self.x_step, 
                                                   self.wavelength,
                                                   elements)
        #formfactors.plot_formfactors()
        self.q = formfactors.get_q()
        self.theta = formfactors.get_theta()
        
        
        global CaF, FCa, CaCa, FF, CeO, OCe, CeCe, FeFe, FeO, OFe, OO
        
        CaF = formfactors.get_formfactor("Ca2+") * formfactors.get_formfactor("F1-") 
        FCa = CaF
        CaCa = formfactors.get_formfactor("Ca2+") * formfactors.get_formfactor("Ca2+")
        FF = formfactors.get_formfactor("F1-") * formfactors.get_formfactor("F1-")
        CeO = formfactors.get_formfactor("Ce4+") * formfactors.get_formfactor("O2-") 
        OCe = CeO
        CeCe = formfactors.get_formfactor("Ce4+") * formfactors.get_formfactor("Ce4+")
        FeFe = formfactors.get_formfactor("Fe3+") * formfactors.get_formfactor("Fe3+")
        FeO = formfactors.get_formfactor("Fe3+") * formfactors.get_formfactor("O2-")
        OFe = FeO
        OO = formfactors.get_formfactor("O2-") * formfactors.get_formfactor("O2-")
        
        print("\nDSE calculator initialized")
        print("\t 2\u03B8-range: %.2f-%.2f" % (self.theta.min()*2, self.theta.max()*2))
        print("\t Q-range: %.2f-%.2f" % (self.q.min(), self.q.max()))
        print("\t wavelength: %.5f \u212B" % (self.wavelength))
    

    def partition_crystal(self):
        """
        Slice the crystal.
        
        Divides the distance and f_ij arrays into parts. The number of the 
        partitions is specified during initialization of the DseCalculator
        with the num_partitions parameter.

        Returns
        -------
        sliced_crystal : ndarray
            Array containing the distances in num_partitions separate arrays.
        sliced_formfactors : chararray
            Array containing the str labels for the formfactor products in 
            num_partitions separate chararrays.

        """
        num_distances = len(self.distances)
        step_size = int(num_distances/self.num_partitions)
        sliced_crystal = []
        sliced_formfactors = []
    
        for i in range(self.num_partitions):
            if i == self.num_partitions-1:
                sliced_crystal.append(self.distances[step_size*i:])
                sliced_formfactors.append(self.f_ij[step_size*i:])
            else:
                sliced_crystal.append(self.distances[step_size*i:step_size*(i+1)])
                sliced_formfactors.append(self.f_ij[step_size*i:step_size*(i+1)])
        return sliced_crystal, sliced_formfactors



    def calculate_Dse(self, distances, f_ij):
        """
        Calculate the Debye scattering equation for non-equal indices.
        
        The str labels of the formfactor products are evaluated with the built-in
        function eval(). f_ij_e now contains f(q) arrayis for every pair of indices.
        f_ij_e has shape (num_pairs, q_steps).
        Distances in angstrom are calculated by multiplying the distances 
        obtained previously with the lattice parameter. This distance array is
        reshaped into a column vector with the numpy method reshape() in order
        to allow broadcasting with f_ij_e array and q array. q*rij_reshaped 
        results in an array with shape (num_pairs, q_steps).
        The Debye scattering equation is evaluated for the upper triangle of indices
        using the numexpr module, afterwards all I_ij(q) arrays are summed up
        to obtain I(q). Numexpr is a fast numerical expression evaluator for NumPy 
        developed by David M. Cooke et al., making use of all cores to do as many
        calculations in parallel as possible. Additionally better performance than
        NumPy is achieved by avoiding allocation of memory for intermediate results
        and making use of a virtual machine written in C.
        For more information see the documentation 
        http://numexpr.readthedocs.io/en/latest/ .

        Parameters
        ----------
        distances : ndarray
            Array containing the pair wise distances.
        f_ij : chararray
            Array containing the formfactor product str labels in same order as
            the distances.

        Returns
        -------
        ndarray
            Scattered intensity of the sum over non-equal indices.

        """
        
        f_ij_e = []
        aa = eval(list((self.atom_numbers.keys()))[0]*2)
        ab = eval(list((self.atom_numbers.keys()))[0]+list((self.atom_numbers.keys()))[1])
        bb = eval(list((self.atom_numbers.keys()))[1]*2)
        for i in f_ij:
            if i == 5:
                f_ij_e.append(ab)
            elif i == 6:
                f_ij_e.append(bb)
            else:
                f_ij_e.append(aa)
        
            #f_ij_e.append(eval(i))
            
        f_ij_e = np.array(f_ij_e)       
        distances *= self.latticepar
        rij_reshaped = distances.reshape(len(distances),1)    
        q = self.q
   
        fijsincqr = ne.evaluate("f_ij_e*sin(q*rij_reshaped)/(q*rij_reshaped)")
        I_sum = fijsincqr.sum(axis=0)

        return 2*I_sum

    def calculate_Dse_total(self):
        """
        Calculate the Debye scattering equation for the whole nanocrystal.
        
        The partitioned distance and formfactor arrays are retrieved from the
        partition_crystal method. The for loop is over the number of partitions
        set during initialization. The intensities obtained from the parts are
        summed up and stored in I_ne (intensity for non-equal indices). The 
        contributions from pairs with equal indices is obtained by summing over
        the kinds of atoms and adding the product of the number of atoms and the
        squared formfactor to I_diagonal. The numbers of atoms are stored in the
        atom_numbers dictionary. The squared formfactors are retrieved from the 
        global variables corresponding to the str label. The string label is 
        obtained by doubling the respective element label, e.g. eval(i*2) with 
        i = 'Fe' --> eval('FeFe') = FeFe = ndarray.

        Returns
        -------
        ndarray
            Total scattered intensity as a function of q.

        """  
        sliced_crystal, sliced_formfactors = self.partition_crystal()
    
        I_ne = 0
        for i in range(self.num_partitions):
            I_ne += self.calculate_Dse(sliced_crystal[i], sliced_formfactors[i])
            print("\t\t %d/%d finished" % (i+1, self.num_partitions))
    
        I_diagonal = 0    
        for i in self.atom_numbers.keys():
            I_diagonal += self.atom_numbers[i]*eval(i*2)
    
        self.intensity = I_diagonal + I_ne
        return self.intensity
       

    def save_to_text(self, filename):
        """
        Save the data.
        
        Saves the data in ASCII columns file, with columns q, 2Theta, I

        Parameters
        ----------
        filename : str
            The filename under which the data will be saved. If no output path
            is given the file will be stored in the package folder.

        Returns
        -------
        None.

        """
        data = np.column_stack((self.q,self.theta*2,self.intensity))
        header = "CIF-file: "+ self.cif_filename
        header += "\nOccupancies: "
        for key in self.occupancies.keys():
            header += key + ": " + str(self.occupancies[key])+ " "
        header += "\nCrystal shape: " + self.shape
        if self.shape == "cube":
             header += "\nEdgelength: %.2f nm" % (self.diameter)
        else:
            header += "\nDiameter: %.2f nm" % (self.diameter)
        
        header += "\nWavelength: %.5f \u212B \nq\t\t\t\t2Theta\t\t\tIntensity\n" % (self.wavelength)
        np.savetxt(self.output_path+filename, data, header=header)
        print("\nFile saved: " + self.output_path+filename)
        
    def perform_calculation(self, filename=None, plot=None):
        """
        Wrapper function to ensure results are saved and plotted or intentionally 
        discarded.

        Parameters
        ----------
        filename : str, optional
            Filename used for saving the data. The default is None.
        plot : bool, optional
            Bool flag for plotting the result. The default is True.

        Returns
        -------
        None.

        """
        print("\n\t starting calculation...")
        I = self.calculate_Dse_total()
        print("\t calculation finished.")
        if plot != None:
            if plot == "Q":
                plotting.plot_results("Q",self.q, self.intensity)
            elif plot == "2Theta":
                plotting.plot_results("2Theta",self.theta*2, self.intensity)
        if filename != None:
            self.save_to_text(filename)
        if filename == None:
            input1 = input("Do you want to save the file? [y/n]\n")
            if input1 == "y":
                input2 = input("Please enter a filename: ")
                self.save_to_text(input2)
            elif input1 == "n":
                return
        

    

