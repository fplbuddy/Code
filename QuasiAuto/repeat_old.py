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
sys.path.insert(1, '/mi/share/scratch/winchester/ForAuto')
import Funcs


# inputs
compname = 'Polaris'
path = os.path.abspath(os.getcwd()) # Getting the path base, might hard code this eventually
paramspath = '/mi/share/scratch/winchester/ForAuto/params.txt'
sPath = '/mi/share/scratch/winchester/ForAuto/SIC'  # path where Shearing ICs are stored THIS NEEDS TO CHANGE
nsPath = '/mi/share/scratch/winchester/ForAuto/NSIC' # path where Non-Shearing ICs are stored THIS NEED TO CHANGE
thresh = 60
cost = 12
CodePath = 'RBC2D_v2.14'
waitTime = 60*60*2 # two hours

            
while True:
    KeepGoing = Funcs.MakeRun(compname,path,paramspath,sPath,nsPath,thresh,cost,CodePath)
    time.sleep(waitTime)
