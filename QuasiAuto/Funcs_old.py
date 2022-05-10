#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 20:28:19 2020

@author: philipwinchester
"""
import numpy as np
import shutil
import math
import time
from datetime import datetime
import os

def RatoRaS(Ra):
    RaChar = str(Ra)
    RaChar = RaChar.replace('.','')
    power =  str(int(np.floor(np.log10(Ra))))
    dp = len(RaChar.rstrip("0"))
    RaStart = RaChar[0]
    if dp == 1:
        return 'Ra_' + RaStart + 'e' + power 
    else:
        RaEnd = RaChar[1:dp]
        return 'Ra_' + RaStart + '_' +  RaEnd +'e' + power 
    
def PrtoPrS(Pr):
    PrChar = str(Pr)
    return 'Pr_' + PrChar

def MakeRun(compname,path,paramspath,sPath,nsPath,thresh,cost,CodePath):
    # Loading in the params file and picking out params
    f=open(paramspath,"r")
    lines=f.readlines()
    f.close()
    # Checking if we have done all we need and then close job
    isFull  =  len(lines[-1]) > 18
    if isFull:
        print('params.txt is full')
        print('')
        return
    load = os.getloadavg()
    load = load[0]
    now = datetime.now()
    current_time = now.strftime("%Y/%m/%d, %H:%M:%S")
    print('The load average was ' + f"{load:.2}" + '% at ' + current_time)
    if load < thresh - cost:
        numRun = math.floor((thresh-load)/cost) # make sure this is not negative
        print('So there is room for ' + str(numRun)  + ' job(s)' )
        # setting up lists we will save things in
        SoNS_list = []
        Pr_list = []
        Ra_list = []
        for i,x in enumerate(lines): # looping round lines in params.txt
            if len(x) < 18: # then we probably have not done it yet
                numRun -= 1 # We have taken one, so remove one from here
                xs = x.split(' ')
                Ra = float(xs[0])
                Pr = float(xs[1])
                ty = xs[2].replace('\n', '')
                SoNS_list.append(ty)
                Ra_list.append(Ra)
                Pr_list.append(Pr)
                lines[i] = x[:-1] + ' ' + compname + '\n'
            if numRun == 0: # If we have got all we need
                break
        # printing the jobs that we are going to submit
        print('they are')
        p = 'Ra = '
        for Ra in Ra_list:
            p = p + "{:.2e}".format(Ra) + ', '
        p = p[0:-2]
        print(p)
        p = 'Pr = '
        for Pr in Pr_list:
            p = p + str(Pr) + ', '
        p = p[0:-2]
        print(p)
        p = 'ty = '
        for ty in SoNS_list:
            p = p + ty + ', '
        p = p[0:-2]
        print(p)

        # updating params.txt
        f = open(paramspath, "w")
        new_file_contents = "".join(lines)
        f.write(new_file_contents)
        f.close()
        for i in range(0,len(Ra_list)):
            Ra = Ra_list[i]
            RaS = RatoRaS(Ra)
            Pr = Pr_list[i]
            PrS = PrtoPrS(Pr)
            SoNS =  SoNS_list[i]
            pathi = path + '/' + SoNS
            # Make S/N path if needed
            if not(os.path.isdir(pathi)): # Checking if path exists, if not make it
                os.mkdir(pathi)
            # Make Pr path if needed
            if not(os.path.isdir(pathi + '/' + PrS)): # Checking if PrS path exists, if not, make it.
                os.mkdir(pathi + '/' + PrS)
            
            pathii = pathi + '/' + PrS + '/' + RaS # This is the final path
            os.mkdir(pathii)
            
            sDirs = ['Checks', 'Fields', 'Spectra']
            for Dir in sDirs:
                os.mkdir(pathii + '/' + Dir) # Adding the dirs where we save the integration
            # Copying things over which i need for the run
            shutil.copy(os.path.join(path,CodePath, 'cleanall.sh'), pathii)
            shutil.copy(os.path.join(path,CodePath, 'input.prm'), pathii)
            shutil.copy(os.path.join(path,CodePath, 'status.prm'), pathii)
            shutil.copy(os.path.join(path,CodePath, 'run.exe'), pathii + '/' + PrS + '_' + RaS + '_' + SoNS + '.exe')
            # copying ICs
            if SoNS == 'S':
                ICPath = sPath
            else:
                ICPath = nsPath 
            src_files = os.listdir(ICPath)
            for file_name in src_files:
                full_file_name = os.path.join(ICPath, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, os.path.join(pathii,'Fields'))
            # Change thing in input file
            inputprm = open(os.path.join(pathii,'input.prm'))
            inputprm_list = inputprm.readlines()
            inputprm.close()
            # Removing start spaces in string
            for i,sentance in enumerate(inputprm_list):
                inputprm_list[i] = sentance.strip()
            # editing the rows
            inputprm_list[0] = pathii + inputprm_list[0] # Fields
            inputprm_list[1] = pathii + inputprm_list[1] # Check
            inputprm_list[2] = pathii + inputprm_list[2] # Spectra
            Rainput = "{:.2e}".format(Ra) # Putting Ra in the format we want
            Rainput.replace('+', '')
            if Rainput[-2] == '0':
                Rainput = Rainput[0:len(Rainput)-3] + Rainput[-1]
            inputprm_list[3] = Rainput + inputprm_list[3][5:]  # Ra 
            Prinput = "{:.2f}".format(Pr)
            inputprm_list[4] = Prinput + inputprm_list[4][4:]  # Pr
            inputprm_list[17] = '0' + inputprm_list[17][1:] # InC
            inputprm_list[19] = '0' + inputprm_list[19][1:] # pert
            inputprm_list[20] = '0' + inputprm_list[20][1:] # norm
            # add back start spaces
            for i,sentance in enumerate(inputprm_list):
                inputprm_list[i] = '  ' + inputprm_list[i] + '\n'
            # save input.prm
            inputprm = open(os.path.join(pathii,'input.prm'), "w")
            new_file_contents = "".join(inputprm_list)
        
            inputprm.write(new_file_contents)
            inputprm.close()
            os.chdir(pathii)
            os.system('mpirun -np 8 ./*.exe &') # running the job
            time.sleep(15) # wait 15 seconds after submitting a job
            os.chdir(path)
    else:
        print('So there is no room to submit job')
    print('')

