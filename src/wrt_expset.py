#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#libralies
import os
import itertools
import numpy as np
import sys
import errno
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
import datetime
import functools
import numpy.random as rd
import os.path
import datetime as dt
import glob
import shutil
import scipy.linalg as spla
from numpy import ma
import random
import re
import calendar
import math

#external python codes
import params as pm
###########################
def assimlation_mode(conflag):
    #  1 - Directly values 
    #  2 - Anomalies
    #  3 - Normalized values
    #  4 - Log converted values
    if conflag==1:
        return "Direct Value"
    if conflag==2:
        return "Anomalies"
    if conflag==3:
        return "Normalized values"
    if conflag==4:
        return "Log converted values"
###########################
def patch_character(patch_size):
    if patch_size==100:
        return "Emperical Local Patch"
    else:
        return "%d"%(patch_size)
###########################
def inflation_para(rho):
    if rho==-1.0:
        return "Adaptive Inflation"
    else:
        return "Inflation Parameter %3.2f"%(rho)
###########################
def calibration(cal):
    if cal=="yes":
        return "calibrated (Xudong et al,. 2021)"
    elif cal=="corrupt":
        return "corrupted (Revel et al,. 2022)"
    else:
        return "not calibrated (Yamazaki et al,. 2011)"
###########################
def corruption(corrupt):
    if corrupt==0:
        return "not corrupted"
    elif corrupt==1:
        return "rivhgt corrupted"
    elif corrupt==2:
        return "rivwth corrupted"
    elif corrupt==3:
        return "rivman corrupted"
    elif corrupt==4:
        return "fldhgt corrupted"
    elif corrupt==5:
        return "all parameters corrupted"
    else:
        return "not corrupted"
###########################
def stat_name(conflag,cal):
    if conflag == 1:
        meanname="None"
        stdname="None"
    if conflag == 2:
        meanname=pm.stat_name(cal)
        stdname="None"
    if conflag == 3:
        meanname=pm.stat_name(cal)
        stdname=pm.stat_name(cal)
    if conflag == 4:
        meanname="None"
        stdname="None"
    return [meanname, stdname]
###########################
def write_text():
    with open("./experimetal_settings.log", "w") as f:
        f.write("# Experimental Settings\n")
        f.write("# "+pm.version()+"\n")
        f.write("# "+pm.CaMa_ver()+"\n")
        f.write("#======================================================\n")
        # Experimental Settings
        f.write("# Experiment Name: "+pm.experiment()+"\n")
        f.write("# Experiment Mode: "+"%d"%(pm.mode())+"\n")
        # Runoff data
        f.write("# Runoff Data: "+pm.input(pm.mode())+", "+pm.runoff_dir()+"\n")
        # Observations
        f.write("# Observations: "+pm.obs_name()+", "+pm.obs_dir()+"\n")
        # Time domain for analysis
        f.write("# Start Date: %04d-%02d-%02d\n"%(pm.starttime()))
        f.write("# End Date: %04d-%02d-%02d\n"%(pm.endtime()))
        # Calibration/Corruption
        f.write("# Model Calibration: "+calibration(pm.calibrate())+"\n")
        # f.write("# Parameter Corruption: "+corruption(pm.corrupt())+"\n")
        # Assimilation Settings
        f.write("# Assimilation Mode: "+assimlation_mode(pm.conflag())+"\n")
        f.write("# Assimilation Domain: \n")
        f.write("# \tWest : %5.2f\n"%(pm.assimW()))
        f.write("# \tSouth: %5.2f\n"%(pm.assimS()))
        f.write("# \tEest : %5.2f\n"%(pm.assimE()))
        f.write("# \tNorth: %5.2f\n"%(pm.assimN()))
        f.write("# Assimilation Settings: \n")
        f.write("# \tPatch Size : "+patch_character(pm.patch_size())+"\n")
        f.write("# \tPatch Id : "+pm.patch_name()+" "+pm.patch_id()+"\n")
        f.write("# \tInflation Method : "+inflation_para(pm.rho())+"\n")
        f.write("# Assimilation Statistics: \n")
        f.write("# \tMean : "+stat_name(pm.conflag(),pm.calibrate())[0]+"\n")
        f.write("# \tStandrad Deviation : "+stat_name(pm.conflag(),pm.calibrate())[1]+"\n")
        f.write("# Created at : "+str(datetime.datetime.now()))
    return 0
###########################
if __name__=="__main__":
    write_text()