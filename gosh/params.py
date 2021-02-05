#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re

########################
#
# parameters list
#
########################

def timestep():
    return 86400 # outer timestep in seconds

def starttime():
    return (2002,1,1) # start date: [year,month,date]

def endtime():
    return (2004,1,1) # end date: [year,month,date]
                      # *note: this date is not included
def spinup_end_year():
    return 2001

def spinup_end_month():
    return 12

def spinup_end_date():
    return 31

def mode():
    return 3
    # parameter to change assimilation mode
    # runoff ensembles will change accordingly.
    # 1: Earth2Obs, 2: ERA20CM, 3: VIC_BC, 4: -25% baised (ELSE_KIM2009/E2O/ERA20CM)

def runname(num):
    if num == 1:
        return "E2O"

    if num == 2:
        return "ERA20CM"

    if num == 3:
        return "VIC_BC"

    if num == 4: #biased runoff experiment
        #return "ELSE_KIM2009"
        return "E2O"
        #return "ERA20CM"

def input(num=mode()):
    if num==1:
        return "E2O"

    if num==2:
        return "ERA20CM"

    if num==3:
        return "VIC_BC"
    # define the runoff data type.

def experiment():
    f=open("./exp.txt","r")
    line=f.readline()
    exp =line.split("\n")[0]
    f.close()
    return exp
    
def ens_mem(mode=mode()):
    if mode == 1:
        return 21
    
    if mode == 2:
        return 20

    if mode == 3:
        return 20

    if mode == 3:
        return 20
    # number of ensemble members

def max_lat():
    return 80. # maximum latitude of assimilation
               # *note: SWOT ovservation is not available beyond 80 degs. this should be less or equal to 80
               ## modified 2018-06-05

def distopen(num):
    if num == 1:
        return 1.0

    if num == 2:
        return 1.0

    if num == 3:
        return 1.0
    #return 0.75 # not needed for ERA20CM
    # corrupted runoff's percentage
    # 0.75 for original Data Assimilation simulation (25% reduced)
    # 1.25 for 25% increased simulation
    # 1.00 for simulation using 1 year before runoff
    # *note: also editing and and re-compile of control_inp at CaMa-Flood is nessessary

def diststd(num):
    if num == 1:
        return 0.1

    if num == 2:
        return 0.25

    if num == 3:
        return 0.25

    #return 1.0 # not needed for ERA20CM
    # noise to make runoff input to scatter ensembles

def assimS():
    return -20
    #return -75
    # data Assimilation's Region (South Edge at latitude)
    # *note: should be larger or equal to -80

def assimN():
    return 5
    #return 75
    # data Assimilation's Region (North Edge at latitude)
    # *note: should be smaller or equal to 80

def assimW():
    return -80
    #return -170
    #return -68.25 # use this for disabling west side of the Amazon basin's observation
    # data Assimilation's Region (West Edge at latitude)
    # *note: should be larger or equal to -170

def assimE():
    return -45
    #return 170
    # data Assimilation's Region (East Edge at latitude)
    # *note: should be smaller or equal to 170

def patch_size():
    #return 0
    return 100
    # the size of the local patch of LETKF(Local ** EnKF)
    # 0: only 1 pixel (the pixel itself) belongs to its local patch
    # 100: empirical local patch

def err_expansion():
    return 1.0
    # variance-covariance expansion
    # works well with 0.04

def rivman_error():
    return 0
    #define the experiment with or without rivman error
    # 0 : with out manning error
    # 1 : with manning error: Manning's n depend on river width
    # 2 : with manning error: Manning's n depend spatial covarience
    # 3 : with manning error: Manning's n randomly distributed subbasin
    # 4 : with manning error: Manning's n randomly distributed
    # 5 : with manning error: Manning's n depend on rivseq
    # 6 : with manning error: Manning's n depend on uparea
    # 7 : with manning error: Manning's n depend on

def run_flag():
    return 0
    # 0 run all simulations
    # 1 run only corrupted and assimilated simulations
    # 2 run only true and assimilated simulations
    # 3 run only assimilated simulation

def true_run(num):
    if num == 1:
        return 3 # ecmwf as true

    if num == 2:
        return 4 # ERA04 as true

    if num == 3: # only one will be used
        return 3

def CaMa_dir():
    return "/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
    #return "/cluster/data6/menaka/CaMa-Flood_v396_20191225"
    #return "/cluster/data6/menaka/CaMa-Flood_v395b_20191030"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders

def mapname():
    return "amz_06min"

def map_dimension():
    fname=CaMa_dir()+"/map/"+mapname()+"/params.txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    #-------
    nx     = int(filter(None, re.split(" ",lines[0]))[0])
    ny     = int(filter(None, re.split(" ",lines[1]))[0])
    gsize  = float(filter(None, re.split(" ",lines[3]))[0])
    return nx,ny,gsize

def DA_dir():
    return "/cluster/data6/menaka/HydroDA"
    # directory of HydroDA
    # where src, dat, sat, out exsits

def patch_dir():
    return "/cluster/data6/menaka/Empirical_LocalPatch/local_patch"
    #return "/cluster/data6/menaka/covariance/local_patch"
    #return "/cluster/data6/menaka/covariance/local_patchMS"
    #return "/cluster/data6/menaka/covariance/local_patch_0.80"

def patch_name():
    return "amz_06min_S14FD"

def patch_id():
    return "0.60"

def spinup_mode():
    return 0
    # 0: do spinup simulation for both (corrupted and true) simulation
    # 1: do spin up only at corrupted simulation
    # 2: do spin up only at true simulation
    # 3: no spinup simulation at all
    ### if initial restart file is ready, spinup simulation is no need

def ovs_err():
    return 0.100000
    # size of SWOT observation error in meters
    # 10cm for 1km^2 water body or 25cm < 1km^2 water area
    # should be at least 0.02
    # hope to be below 0.10

def thersold():
    return 0.60
    # thersold to define the local patch

def initial_infl():
    return 1.08
    # initial inflation parameter

def rho():
    return -1.0
    # -1.0 : adaptive inflation will be used Myoshi et al (2011)
    # positive : fixed inflation parameter will be used

def sigma_b():
    return 0.0400000
    # bacground variance of inflation for adaptive inflation Myoshi et al (2011)

def MKLdir():
    return "/opt/intel/compilers_and_libraries_2016.3.170/mac/mkl"
    # directory of Intel MKL files
    # Intel MKL is needed for doing data assimilation
    # Please Download and Instal it to your System before running
    # for more information --> https://software.intel.com/en-us/qualify-for-free-software/academicresearcher

def output_er():
    return 0
    # setting for saving or deleting intermediate files
    # 0 for saving & 1 for deleting
    # those files may be more than 400GB, so erasing is recommended if not necessary

def HydroWeb_dir():
    return "/cluster/data6/menaka/HydroWeb"

def make_log():
    return 1
    # setting for making log files
    # 1 is for making and 0 is for not making

def para_nums():
    return 10
    # setting number of parallels to run CaMa-Flood Model
    # defualt is 6, but may change depending on your system

def slack_notification():
    return 0
    # setting for validating slack notification
    # 1 for valid and 0 for invalid
    # 0 is a default if you are not familiar with slack
    # if you turn it to 1, you need to edit sendslack.py
    # for more information refer https://api.slack.com/incoming-webhooks

def ens_at_non():
    return 1
    # * At Recent version, ensemble generating random number is constant for full simulation.
    # (For example, when ensemble 001 is corrupted with -0.1 at day 1, ensemble 001 will be always corrupted with 0.1 for full simulation period.)
    # Previously, ensemble mean was used as an assimilated value for non-observed location.
    # In this version, this treatment has changed and non-observed location is given with an ensemble value.
    # To enable this new feature, set the return of params.py method “ens_at_non()”, “1”(DEFAULT SETTING).
    # If you don’t want to use this, set it to “0”.

# functions for corrupting manning coeffcient ###################
# this is for corrupting manning coefficient at Corrupted Simulation
# manning coefficient will be corrupted with random numbers generated from following functions
# the random number is generated for each ensemble member
# random number is made by gaussian noise of average = corruptman_base(), stddev = corruptman_std()
def corruptman_base():
    return 0.03

def corruptman_std():
    return 0.015

def rivman_base():
    return 0.03

def rivman_min():
    return 0.025

def rivman_max():
    return 0.035

def corruptele_base():
    return 0.5 # not needed for ERA20CM

def corruptele_std():
    return 0.25 # not needed for ERA20CM

def non_hgt():
    return 7.0 # not needed for ERA20CM
    # nominal water height

def cpu_nums():
    f=open("./ncpus.txt","r")
    line=f.readline()
    ncpus =int(line.split("\n")[0])
    f.close()
    return ncpus/para_nums()
    # number of cpus used

def version():
    return "v5.0.0 CaMa-Flood and HydroWeb"#, Bathymetry"
    # version for WSE assimilation / observation localization
    # CaMa-Flood v396a used