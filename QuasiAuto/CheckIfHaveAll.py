#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 17:12:30 2021

@author: philipwinchester
check if have all
"""
import os
import FuncsQuasi as Funcs

SamsungPath = '/Volumes/Samsung_T5/Quasi3D'

# Reads in the file with the runs which we want to move over
path = os.getcwd()
f=open(path + '/paramsQuasi.txt',"r")
lines=f.readlines()
f.close()

for i,x in enumerate(lines): # looping round lines paramsQuasi_want.txt
    xs = x.split(' ')
    Ra = float(xs[0])
    Pr = float(xs[1])
    SoN = xs[2]
    Ly = float(xs[3].replace('\n', ''))
    # build path
    RaS = Funcs.normaltoS(Ra,'Ra')
    PrS = Funcs.normaltoS(Pr,'Pr')
    LyS = Funcs.normaltoS(Ly,'Ly')
    pathadd = '/IC_'+SoN + '/N_256x256' + '/' + PrS + '/' + RaS + '/' + LyS # List of path chain we want to add
    thingsToAdd = ['/IC_'+SoN, '/N_256x256', '/' + PrS, '/' + RaS, '/' + LyS]
    # Check if we have data in the folder from cluster
    if not(os.path.exists(SamsungPath + pathadd)):
        print(pathadd)
    
            
            
# updating paramsQuasi_want.txt
#f = open(path + '/paramsQuasi_want.txt', "w")
#new_file_contents = "".join(lines)
#f.write(new_file_contents)
#f.close()

