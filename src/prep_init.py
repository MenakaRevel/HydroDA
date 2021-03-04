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
# ################################################
#
# make folders
# Prepare input data sets [runoff]
# initial inflation parameter rho for assimilation
#
# ################################################
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
###########################
def slink(src,dst):
    try:
        os.symlink(src,dst)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src,dst)
        else:
            raise
###########################
def initial(): #used
    # program for initialization

    # creating output folders
    mkdir("CaMa_out")
    mkdir("CaMa_in")
    mkdir("CaMa_in/restart")
    mkdir("CaMa_in/restart/assim")
    mkdir("CaMa_in/restart/open")
    mkdir("CaMa_in/restart/true")
    mkdir("assim_out")
    mkdir("assim_out/xa_m")
    mkdir("assim_out/xa_m/assim")
    mkdir("assim_out/xa_m/open")
    mkdir("assim_out/xa_m/true")

    mkdir("assim_out/ens_xa")         # NEW v.1.1.0
    mkdir("assim_out/ens_xa/assim")   # NEW v.1.1.0
    mkdir("assim_out/ens_xa/open")    # NEW v.1.1.0

    mkdir("assim_out/nonassim")
    mkdir("assim_out/nonassim/open")
    mkdir("assim_out/nonassim/assim")

    mkdir("assim_out/rest_true")
    mkdir("assim_out/rivout")
    mkdir("assim_out/rivout/open")
    mkdir("assim_out/rivout/assim")
    mkdir("assim_out/rivout/true")
    mkdir("assim_out/outflw")
    mkdir("assim_out/outflw/open")
    mkdir("assim_out/outflw/assim")
    mkdir("assim_out/outflw/true")
    mkdir("assim_out/fldout")
    mkdir("assim_out/fldout/open")
    mkdir("assim_out/fldout/assim")
    mkdir("assim_out/fldout/true")
    mkdir("assim_out/flddph")
    mkdir("assim_out/flddph/open")
    mkdir("assim_out/flddph/assim")
    mkdir("assim_out/flddph/true")
    mkdir("err_out")
    mkdir("assim_out/fldarea/")
    mkdir("assim_out/fldarea/open")
    mkdir("assim_out/fldarea/assim")
    mkdir("assim_out/fldarea/true")
    
    # inflation parameter
    mkdir("inflation")

    # observations folder
    # inside the assim_out folder
    mkdir("assim_out/obs")

    os.system("touch assim_out/__init__.py")
    mkdir("logout")

    #mkdir("assim_out/rivhgt")
    #mkdir("assim_out/rivdph")
    #mkdir("assim_out/rivhgt/assim")
    #mkdir("assim_out/rivdph/assim")
    #mkdir("assim_out/rivdph/open")
    #mkdir("assim_out/rivdph/true")

    # link input files
    #os.system("ln -s "+pm.CaMa_dir()+"/inp/ELSE_GPCC ./CaMa_in")
    #slink(pm.CaMa_dir()+"/inp/ELSE_GPCC","./CaMa_in/ELSE_GPCC")

    #mkdir("./CaMa_in/ELSE_GPCC/mean_month")
    #mkdir("./CaMa_in/ELSE_KIM2009/mean_month")

    # make ./CaMa_in/ELSe_GPCC/* runoff files
    #mkdir("./CaMa_in/ELSE_GPCC/Roff_TRUE")
    #mkdir("./CaMa_in/ELSE_GPCC/Roff_CORR")
    #mkdir("./CaMa_in/ELSE_KIM2009/Roff_TRUE")
    #mkdir("./CaMa_in/ELSE_KIM2009/Roff_CORR")


    return 0
###########################
def make_initial_infl():
    nx,ny,gsize=pm.map_dimension()
    parm_infl=np.ones([ny,nx],np.float32)*pm.initial_infl()
    start_year,start_month,start_date=pm.starttime() # Start year month date
    yyyy='%04d' % (start_year)
    mm='%02d' % (start_month)
    dd='%02d' % (start_date)
    parm_infl.tofile(pm.DA_dir()+"/out/"+pm.experiment()+"/inflation/parm_infl"+yyyy+mm+dd+".bin")
    return 0
###########################
def make_anomaly_data(mode=pm.mode()):
    # make directory for mean sfcelv
    mkdir(pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/")
    # copy the anomaly files
    if mode == 1:
        # for mean
        iname=pm.DA_dir()+"/dat/mean_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname=pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/mean_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
        # for std
        iname=pm.DA_dir()+"/dat/std_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname=pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/std_sfcelv.bin"
        os.system("cp "+iname+" "+oname)

    if mode == 2:
        # for mean
        iname=pm.DA_dir()+"/dat/mean_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname=pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/mean_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
        # for std
        iname=pm.DA_dir()+"/dat/std_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname=pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/std_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
    
    if mode == 3:
        # for mean
        iname=pm.DA_dir()+"/dat/mean_sfcelv_VIC_BC_"+pm.mapname()+"_1979-2013.bin"
        oname=pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/mean_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
        # for std
        iname=pm.DA_dir()+"/dat/std_sfcelv_VIC_BC_"+pm.mapname()+"_1979-2013.bin"
        oname=pm.DA_dir()+"/out/"+pm.experimet_name()+"/assim_out/mean_sfcelv/std_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
    return 0
###########################
def save_statistic():
    # copy mean and std of simulated WSE
    # for anomaly and normalized assimilations
    mkdir("./assim_out/mean_sfcelv/")
    if pm.input()=="E2O":
        os.system("cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_2000-2014.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin")
        os.system("cp -r "+pm.DA_dir()+"/dat/std_sfcelv_"+pm.input()+"_"+pm.mapname()+"_2000-2014.bin ./assim_out/mean_sfcelv/std_sfcelv.bin")
        print "cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_2000-2014.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin"
    if pm.input()=="VIC_BC":
        os.system("cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin")
        os.system("cp -r "+pm.DA_dir()+"/dat/std_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin ./assim_out/mean_sfcelv/std_sfcelv.bin")
        print "cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin"
    return 0
###########################
# # make necessary directories
# print "initial"
# initial()

# # prepare runoff ensembles
# print "prepare input"
# prepare_input()

# # initial inflation parameter rho for assimilation
# print "make intial inflation"
# make_initial_infl()

# # preapre the mean and std for anomaly/normalized assimilation
# print "save statistics"
# save_statistic()