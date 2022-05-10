#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:30:57 2021

@author: philipwinchester
"""

#from apscheduler.schedulers.blocking import BlockingScheduler
import os
os.path.realpath(__file__) # This moves to where the file is
import time
import sys
#sys.path.insert(1, '/mi/share/scratch/winchester/ForAuto')
sys.path.insert(1, '/Users/philipwinchester/Dropbox/Codes/QuasiAuto')
import FuncsQuasi

# inputs
compname = 'Polaris'
path = os.path.abspath(os.getcwd()) # Getting the path base, might hard code this eventually
#paramspath = '/mi/share/scratch/winchester/ForAuto/paramsQuasi.txt'
paramspath = '/Users/philipwinchester/Dropbox/Codes/QuasiAuto/paramsQuasi.txt'
#ICPath = '/home/winchester/Documents/ICs/Quasi3D/'  # path where Shearing ICs are stored THIS NEEDS TO CHANGE
ICPath = '/Users/philipwinchester/Dropbox/Codes/QuasiAuto/RandomFolder/' 
thresh = 95
#CodePath = 'Q3D_v3'
CodePath = 'test_code'
waitTime = 60*60*2 # two hours


while True:
    KeepGoing = FuncsQuasi.MakeRun(compname,path,paramspath,ICPath,thresh,CodePath)
    time.sleep(waitTime)
