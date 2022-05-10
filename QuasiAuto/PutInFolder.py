#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:22:39 2021

@author: philipwinchester
Some code that takes data from cluster and puts in Samsung
"""


import os
import FuncsQuasi as Funcs
from distutils.dir_util import copy_tree

ClusterPath = '/Users/philipwinchester/Desktop/hello'
SamsungPath = '/Volumes/Samsung_T5/Quasi3D'

# Reads in the file with the runs which we want to move over
path = os.getcwd()
f=open(path + '/paramsQuasi.txt',"r")
lines=f.readlines()
f.close()
add = 0
alreadygot = 0
for i,x in enumerate(lines): # looping round lines paramsQuasi_want.txt
    xs = x.split(' ')
    Ra = float(xs[0])
    Pr = float(xs[1])
    SoN = xs[2]
    Ly = float(xs[3])
    Nx = xs[4]
    Ny = xs[5]
    # build path
    RaS = Funcs.normaltoS(Ra,'Ra')
    PrS = Funcs.normaltoS(Pr,'Pr')
    LyS = Funcs.normaltoS(Ly,'Ly')
    res = f'/N_{Nx}x{Ny}'
    pathadd = '/IC_'+SoN + res + '/' + PrS + '/' + RaS + '/' + LyS # List of path chain we want to add
    thingsToAdd = ['/IC_'+SoN, res, '/' + PrS, '/' + RaS, '/' + LyS]
    
    # Check if we have data in the folder from cluster
    if os.path.exists(ClusterPath + pathadd):
        # Now chech if we have in samsung
        if os.path.exists(SamsungPath + pathadd):
            print('We already have ' + pathadd + ' in Samsung' )
            alreadygot += 1
        else:
            print(f'new: {pathadd}')
            add += 1
            # Add the folders we need
            pathi = SamsungPath
            for i in range(0,len(thingsToAdd)):
                pathi += thingsToAdd[i]
                if not(os.path.isdir(pathi)): # Checking if path exists, if not make it
                    os.mkdir(pathi)
            # At the end of this for loop, pathi is where we put the work
            # Now move over the data
            copy_tree(ClusterPath + pathadd, pathi)
            # And now we add an ok to end of paramsQuasi_want.txt line
            #lines[i] = x[:-1] + ' ok' + '\n'
        
add = str(add)
alreadygot = str(alreadygot)
print(f'We have added: {add}')
print(f'We already had: {alreadygot}')
            
            
# updating paramsQuasi_want.txt
#f = open(path + '/paramsQuasi_want.txt', "w")
#new_file_contents = "".join(lines)
#f.write(new_file_contents)
#f.close()
    

    
    
    

