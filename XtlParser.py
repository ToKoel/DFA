#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:45:51 2020

@author: tobiaskohler
"""

class XtlParser():
    
    def __init__(self, filename =""):
        self.filename = filename
        self.elementlist = list()
        self.coordinatelist = list()
        self.latticepar = 0
  
        with open(self.filename) as fp: 
            counter = 1
            for line in fp:
                if line.split()[0] == "EOF":
                    break
                if counter == 3:
                    self.latticepar = float(line.split()[0])
                if counter >= 8:
                    l = line.split()
                    self.elementlist.append(l[0])
                    self.coordinatelist.append([float(l[1]),float(l[2]),float(l[3])])
                counter +=1
    
                
    