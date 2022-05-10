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
import psutil
import multiprocessing

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


def normaltoS(val,s):
    RaChar = str(val)
    RaChar = RaChar.replace('.','')
    power =  str(int(np.floor(np.log10(val))))
    power = power.replace('-','_')
    while RaChar[0] == "0": # Remove zeros at the start
        RaChar = RaChar[1:len(RaChar)]
    while RaChar[-1] == "0": # Remove zeros at end
        RaChar = RaChar[0:len(RaChar)-1]
    dp = len(RaChar)
    RaStart = RaChar[0]
    if dp == 1:
        return s + '_' + RaStart + 'e' + power
    else:
        RaEnd = RaChar[1:dp]
        return s + '_' + RaStart + '_' +  RaEnd +'e' + power


def StoNormal(val):
    eloc = val.find('e')
    base = val[val.find("_")+1:eloc]
    base = base.replace('_','.')
    power = val[eloc+1:len(val)]
    power = power.replace('_','-')
    return float(base)*10**int(power)


def DNSformat(val):
    Rainput = "{:.2e}".format(val) # Putting Ra in the format we want
    Rainput = Rainput.replace('+', '')
    if Rainput[-2] == '0':
        Rainput = Rainput[0:len(Rainput)-2] + Rainput[-1]
    return Rainput

def ICpathMaker(path,SoN,PrS,RaS,Nx,Nz):
    if SoN == "S":
        path = path + 'N_' + Nx + 'x' + Nz + '/' + PrS + '/Shearing/'
    else:
        path =  path + 'N_' + Nx + 'x' + Nz + '/' + PrS + '/NonShearing/'
    # Find the Ra that we have
    RaS_list = os.listdir(path)
    Ra_list = [] # Will add Ra to this
    for i in range(len(RaS_list)):
        if "Ra_" in RaS_list[i]:
            RaSinst = RaS_list[i]
            Ra = StoNormal(RaSinst)
            Ra_list.append(Ra)
    Ra = StoNormal(RaS)
    log_absolute_difference_function = lambda list_value : max(list_value/Ra, (list_value/Ra)**-1) # Distance function
    RaC = min(Ra_list, key=log_absolute_difference_function)
    RaSC = normaltoS(RaC,"Ra")
    path = path + RaSC
    return path

def GetNumberRuns():
    count  = 0
    for proc in psutil.process_iter():
        try:
            processName = proc.name()
            if "Pr_" in processName:
                count += 1
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return int(count)

def EmptyCores():
    cpu = psutil.cpu_percent(10,True)
    empt = sum(i < 10 for i in cpu)
    return empt

def GetCPU():
    check = None
    for proc in psutil.process_iter():
        try:
            processName = proc.name()
            if "Pr_" in processName:
                check = proc.cpu_percent(interval=1)
                if check < 90:
                    return check, False # There is a process which has too low CPU
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return check,True

def MakeRun(compname,path,paramspath,ICPathbase,thresh,CodePath):
    # Loading in the params file and picking out params
    f=open(paramspath,"r")
    lines=f.readlines()
    f.close()
    # Checking if we have done all we need and then close job
    isFull  =  len(lines[-1]) > 30
    if isFull:
        print('params.txt is full')
        print('')
        return
    # Pring how mant we have running
    running = GetNumberRuns()
    print('Currently there are ' + str(running) + ' cores running')
    MaxCores = math.floor(multiprocessing.cpu_count()/4)
    CPC =  math.ceil(100/multiprocessing.cpu_count())
    # Now check total CPU load
    load = os.getloadavg()
    load = load[0]
    now = datetime.now()
    current_time = now.strftime("%Y/%m/%d, %H:%M:%S")
    print('The load average was ' + f"{load:.2}" + '% at ' + current_time)
    # Picking our first job
    for i,x in enumerate(lines): # looping round lines in params.txt
        if len(x) < 30: # then we probably have not done it yet
            xs = x.split(' ')
            # Pikcing up the parameters
            Ra = float(xs[0])
            Pr = float(xs[1])
            SoN = xs[2]
            Ly = float(xs[3].replace('\n', ''))
            Nx = xs[4]
            Nz = xs[5]
            nmpi = xs[6].replace('\n', '')
            break # break out of for loop
    # Check if we can submit
    while load + int(nmpi)*CPC < thresh and running + int(nmpi) <= MaxCores:
        # Adding to stuff we measure against
        running += int(nmpi)
        load += int(nmpi)*CPC
        # upbdate params file
        lines[i] = x[:-1] + ' ' + compname + '\n'
        f = open(paramspath, "w")
        new_file_contents = "".join(lines)
        f.write(new_file_contents)
        f.close()
        # print some info about run
        print('Set up (Ra, Pr, Ly) = ' + "{:.2e}".format(Ra) + ', ' + str(Pr) + ', ' + str(Ly))
        # Setting up the run
        RaS = normaltoS(Ra,'Ra')
        PrS = normaltoS(Pr,'Pr')
        LyS = normaltoS(Ly,'Ly')
        thingsToAdd = ['/IC_'+SoN, f'/N_{Nx}x{Nz}', '/' + PrS, '/' + RaS, '/' + LyS] # List of path chain we want to add
        pathi = path
        for i in range(0,len(thingsToAdd)):
            pathi += thingsToAdd[i]
            if not(os.path.isdir(pathi)): # Checking if path exists, if not make it
                os.mkdir(pathi)
        # At the end of this for loop, pathi is where we put the work
        sDirs = ['Checks', 'Fields', 'Spectra']
        for Dir in sDirs:
            os.mkdir(pathi + '/' + Dir) # Adding the dirs where we save the integration
        # Copying things over which i need for the run
        CodePathi = CodePath + '_' + Nx + 'x' + Nz # Putting resolution in CodePath
        shutil.copy(os.path.join(path,CodePathi, 'cleanall.sh'), pathi)
        shutil.copy(os.path.join(path,CodePathi, 'input.prm'), pathi)
        shutil.copy(os.path.join(path,CodePathi, 'status.prm'), pathi)
        shutil.copy(os.path.join(path,CodePathi, 'run.exe'), pathi + '/' + PrS + '_' + RaS + '_' + LyS + '_' + SoN + '.exe')
        # copying ICs
        ICPath = ICpathMaker(ICPathbase,SoN,PrS,RaS,Nx,Nz)
        src_files = os.listdir(ICPath)
        for file_name in src_files:
            full_file_name = os.path.join(ICPath, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.path.join(pathi,'Fields'))
        # Change thing in input file
        inputprm = open(os.path.join(pathi,'input.prm'))
        inputprm_list = inputprm.readlines()
        inputprm.close()
        # Removing start spaces in string
        for i,sentance in enumerate(inputprm_list):
            inputprm_list[i] = sentance.strip()
        # editing the rows
        inputprm_list[0] = pathi + inputprm_list[0] # Fields
        inputprm_list[1] = pathi + inputprm_list[1] # Check
        inputprm_list[2] = pathi + inputprm_list[2] # Spectra
        Rainput = DNSformat(Ra)
        Prinput = DNSformat(Pr)
        Lyinput = DNSformat(Ly)
        inputprm_list[3] = Rainput + inputprm_list[3][5:]  # Ra
        inputprm_list[4] = Prinput + inputprm_list[4][5:]  # Pr
        inputprm_list[5] = Lyinput + inputprm_list[5][5:]  # Ly
        # add back start spaces
        for i,sentance in enumerate(inputprm_list):
            inputprm_list[i] = '  ' + inputprm_list[i] + '\n'
        # save input.prm
        inputprm = open(os.path.join(pathi,'input.prm'), "w")
        new_file_contents = "".join(inputprm_list)

        inputprm.write(new_file_contents)
        inputprm.close()
        os.chdir(pathi)
        os.system('mpirun -np ' + nmpi  + ' ./*.exe &') # running the job
        time.sleep(15) # wait 15 seconds after submitting a job
        os.chdir(path)
        # pikcing up new run
        f=open(paramspath,"r")
        lines=f.readlines()
        f.close()
        # Checking if we have done all we need and then close job
        isFull  =  len(lines[-1]) > 30
        if isFull:
            print('params.txt is full')
            print('')
            return  
        for i,x in enumerate(lines): # looping round lines in params.txt
            if len(x) < 30: # then we probably have not done it yet
                xs = x.split(' ')
                # Pikcing up the parameters
                Ra = float(xs[0])
                Pr = float(xs[1])
                SoN = xs[2]
                Ly = float(xs[3].replace('\n', ''))
                Nx = xs[4]
                Nz = xs[5]
                nmpi = xs[6].replace('\n', '')
                break # break out of for loop
    print('')
    
"""def MakeFolders(compname,path,paramspath,ICPathbase,CodePath,numRun):
    # Loading in the params file and picking out params
    f=open(paramspath,"r")
    lines=f.readlines()
    f.close()
    # setting up lists we will save things in
    SoN_list = []
    Pr_list = []
    Ra_list = []
    Ly_list = []
    for i,x in enumerate(lines): # looping round lines in params.txt
        if len(x) < 20: # then we probably have not done it yet
            numRun -= 1 # We have taken one, so remove one from here
            xs = x.split(' ')
            Ra = float(xs[0])
            Pr = float(xs[1])
            ty = xs[2]
            Ly = float(xs[3].replace('\n', ''))
            SoN_list.append(ty)
            Ra_list.append(Ra)
            Pr_list.append(Pr)
            Ly_list.append(Ly)
            lines[i] = x[:-1] + ' ' + compname + '\n'
        if numRun == 0: # If we have got all we need
            break

    # updating params.txt
    f = open(paramspath, "w")
    new_file_contents = "".join(lines)
    f.write(new_file_contents)
    f.close()
    # Now looping round the jobs
    print('The jobs are in the folders')
    for i in range(0,len(Ra_list)):
        Ra = Ra_list[i]
        RaS = normaltoS(Ra,'Ra')
        Pr = Pr_list[i]
        PrS = normaltoS(Pr,'Pr')
        SoN =  SoN_list[i]
        Ly = Ly_list[i]
        LyS = normaltoS(Ly,'Ly')
        thingsToAdd = ['/IC_'+SoN, '/N_256x256', '/' + PrS, '/' + RaS, '/' + LyS] # List of path chain we want to add
        pathi = path
        for i in range(0,len(thingsToAdd)):
            pathi += thingsToAdd[i]
            if not(os.path.isdir(pathi)): # Checking if path exists, if not make it
                os.mkdir(pathi)
        # At the end of this for loop, pathi is where we put the work

        sDirs = ['Checks', 'Fields', 'Spectra']
        for Dir in sDirs:
            os.mkdir(pathi + '/' + Dir) # Adding the dirs where we save the integration
        # Copying things over which i need for the run
        shutil.copy(os.path.join(path,CodePath, 'cleanall.sh'), pathi)
        shutil.copy(os.path.join(path,CodePath, 'input.prm'), pathi)
        shutil.copy(os.path.join(path,CodePath, 'status.prm'), pathi)
        shutil.copy(os.path.join(path,CodePath, 'run.exe'), pathi + '/' + PrS + '_' + RaS + '_' + LyS + '_' + SoN + '.exe')
        # copying ICs
        ICPath = ICpathMaker(ICPathbase,SoN,PrS,RaS)
        src_files = os.listdir(ICPath)
        for file_name in src_files:
            full_file_name = os.path.join(ICPath, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.path.join(pathi,'Fields'))
        # Change thing in input file
        inputprm = open(os.path.join(pathi,'input.prm'))
        inputprm_list = inputprm.readlines()
        inputprm.close()
        # Removing start spaces in string
        for i,sentance in enumerate(inputprm_list):
            inputprm_list[i] = sentance.strip()
        # editing the rows
        inputprm_list[0] = pathi + inputprm_list[0] # Fields
        inputprm_list[1] = pathi + inputprm_list[1] # Check
        inputprm_list[2] = pathi + inputprm_list[2] # Spectra
        Rainput = DNSformat(Ra)
        Prinput = DNSformat(Pr)
        Lyinput = DNSformat(Ly)
        inputprm_list[3] = Rainput + inputprm_list[3][5:]  # Ra
        inputprm_list[4] = Prinput + inputprm_list[4][5:]  # Pr
        inputprm_list[5] = Lyinput + inputprm_list[5][5:]  # Ly
        # add back start spaces
        for i,sentance in enumerate(inputprm_list):
            inputprm_list[i] = '  ' + inputprm_list[i] + '\n'
        # save input.prm
        inputprm = open(os.path.join(pathi,'input.prm'), "w")
        new_file_contents = "".join(inputprm_list)
        inputprm.write(new_file_contents)
        inputprm.close()
        print(pathi)
    print('')"""
