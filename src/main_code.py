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
#import src.letkf_lib as lb
#import src.calc_one_daybef as odb


########################
#
# main program
#
########################

#main code for LETKF
############
## main control function
############
def main_act():

    print (pm.experiment())
    print (pm.version())
    print (pm.runname(pm.mode()))

    # Set Time
    print ("Set Time")
    timestep=pm.timestep() #time step for assimilation
    start_year,start_month,start_date=pm.starttime() # Start year month date
    end_year,end_month,end_date=pm.endtime() # End year month date

    # Spin-up Simulation
    print ("spin up simulation")
    spin_up()

    # make initial restart
    print ("make intial restart")
    make_initial_restart()
    #make_initial_restart_one()

    start_dt=datetime.date(start_year,start_month,start_date)
    #start_dt=datetime.date(2008,8,1)#5,29)
    end_dt=datetime.date(end_year,end_month,end_date)
    days_count=(end_dt-start_dt).days
    #days_count=(end_dt-datetime.date(1991,1,1)).days 

    # run daily simulations
    for day in np.arange(days_count):
        running_dt=start_dt+datetime.timedelta(days=day)
        yyyy='%04d' % (running_dt.year)
        mm='%02d' % (running_dt.month)
        dd='%02d' % (running_dt.day)

        # if dd=="01" and pm.slack_notification()==1:
        #     os.system("source src/sendslack.sh python_notification DA in progress "+yyyy+" "+mm+" "+dd)

        one_day_loop(yyyy,mm,dd,day)

    # clean all intermediate files
    if pm.output_er()==1:
        os.system("rm -Rf ./CaMa_out/"+yyyy+"*")
        os.system("rm -Rf ./CaMa_in/"+pm.input()+"/Roff_CORR/Roff__"+yyyy+"*")

################
## single loop program
################
def one_day_loop(yyyy,mm,dd,day):
    print ("================================ start loop of "+yyyy+" "+mm+" "+dd+" ========================================")
    '''
    # Corrupted Simulation (Open Loop) ###################################################
    if pm.run_flag() == 0 or pm.run_flag() == 1:
        ODM_inputlist=[]

        # set for ensemble simulations
        ens_num=1
        for ens_num in np.arange(1,pm.ens_mem()+1):
            ODM_inputlist.append([yyyy,mm,dd,'%03d'%ens_num,"open"])

        # Run CaMa-Flood Model (ensemble simulations)
        p=Pool(pm.para_nums())
        p.map(one_day_sim,ODM_inputlist)
        p.terminate()

        # copy ensemble corrupted simulation forecasts from CaMa_out (xa)
        cprestart=[]
        for ens_num in np.arange(1,pm.ens_mem()+1):
            cprestart.append([yyyy,mm,dd,ens_num])
        p=Pool(pm.para_nums())
        p.map(copy_corrupted_sfcelv,cprestart)
        p.terminate()

        #copy_corrupted_sfcelv(yyyy,mm,dd,num) # changed to copy sfcelv
        # copy restart files/ no need for recalculation - MODIFIED @ Menaka
        p=Pool(pm.para_nums())
        p.map(copy_corrupted_restart,cprestart)
        p.terminate()

        # make restart MODIFIED v.1.1.0
        #mkrestart=[]
        #for ens_num in np.arange(1,pm.ens_mem()+1):
        #    mkrestart.append([yyyy,mm,dd,"open",'%03d'%ens_num])

        ## Modify the restart file
        #p=Pool(pm.para_nums())
        #p.map(make_restart,mkrestart)
        #p.terminate()
    '''
    # Assimilated Simulation #############################################################
    ODM_inputlist=[]

    # set for ensemble simulations
    for ens_num in np.arange(1,pm.ens_mem()+1):
        ODM_inputlist.append([yyyy,mm,dd,'%03d'%ens_num,"assim"])

    # Run CaMa-Flood Model (ensemble simulations)
    p=Pool(pm.para_nums())
    p.map(one_day_sim,ODM_inputlist)
    p.terminate()

    # make forecasted value for assimilated simulation
    # do assimilation (LETKF)
    data_assim(yyyy,mm,dd)
    #direct_insert(yyyy,mm,dd,day)

    # make restart MODIFIED v.1.1.0
    mkrestart=[]
    for ens_num in np.arange(1,pm.ens_mem()+1):
        mkrestart.append([yyyy,mm,dd,"assim",'%03d'%ens_num])

    ## ***** no need for restart file gernation
    # # Modify the restart file
    # p=Pool(pm.para_nums())
    # p.map(make_restart,mkrestart)
    # p.terminate()

    # store river variable files
    store_out(yyyy,mm,dd)

#    # make rivout @menaka
#    mkrivout=[]
#    for ens_num in np.arange(1,pm.ens_mem()+1):
#        mkrivout.append([yyyy,mm,dd,"assim",'%03d'%ens_num])
#    p=Pool(pm.para_nums())
#    p.map(make_rivout,mkrivout)
#    p.terminate()

    # clean files
    if pm.output_er()==1:
        bef_dt=datetime.date(int(yyyy),int(mm),int(dd))-datetime.timedelta(days=1)
        bef_yyyy='%04d' %bef_dt.year
        bef_mm='%02d' %bef_dt.month
        bef_dd='%02d' %bef_dt.day
        os.system("rm -Rf ./CaMa_out/"+bef_yyyy+bef_mm+bef_dd+"*")


#######################################################################################


############
## main program functions
############
def spin_up(): #used
    # run spin up simulation
    # 1 year spin up for calculating initial value
    # one simulation for true
    # ensmble simulation for open

    dir2=pm.CaMa_dir()
    cpunums = pm.cpu_nums()
    yyyy = "%04d"%(pm.spinup_end_year())
    print (pm.spinup_mode())
    if pm.spinup_mode()==3:
        return 0

    inputlist=[]

#    if pm.spinup_mode()==0 or pm.spinup_mode()==2:
#        inputlist.append([yyyy,"true",'000'])
#        #spinup_loop(inputlist)

    if pm.spinup_mode()==0 or pm.spinup_mode()==1:
        for ens_num in np.arange(1,pm.ens_mem()+1):
            inputlist.append([yyyy,"open",'%03d'%ens_num])

    # Run spinup simulations
    p=Pool(pm.para_nums())
    p.map(spinup_loop,inputlist)
    p.terminate()

    print ("======================= end spinup ==========================")

    return 0
###########################
def spinup_loop(inputlist):
    # Run spinup simulation
    yyyy=inputlist[0]
    loop=inputlist[1]
    ens_num=inputlist[2]
    dir2=pm.CaMa_dir()
    cpunums=pm.cpu_nums()
    mode=pm.mode()
    run_name=pm.runname(mode)
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    mapname=pm.mapname()
    cal=pm.calibrate()
    print  ("%s for %03d"%(loop,int(ens_num)))
    os.system("source "+pm.DA_dir()+"/src/spin_up.sh "+str(yyyy)+" "+str(loop)+" "+ens_num+" "
    +dir2+" "+str(cpunums)+" "+str(run_name)+" "+str(exp_dir)+" "+str(mapname)+" "+str(cal))
    return 0
###########################
def one_day_sim(inputlist):
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    ens_num=inputlist[3]
    looptype=inputlist[4]
    mode=pm.mode()
    run_name=pm.runname(mode)

    # program for running one day model

    bef_dt=datetime.date(int(yyyy),int(mm),int(dd))-datetime.timedelta(days=1)
    bef_yyyy='%04d' %bef_dt.year
    bef_mm='%02d' %bef_dt.month
    bef_dd='%02d' %bef_dt.day

    print ("oneday loop for",yyyy,mm,dd,ens_num,looptype)
    dir2=pm.CaMa_dir()
    #if looptype=="true":
    #    distopen="1.0"
    #else:
    #    distopen=str(pm.distopen())

    print (yyyy+" "+mm+" "+dd+" "+ens_num+" "+dir2+" "+looptype)
    cpunums = pm.cpu_nums()
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    mapname=pm.mapname()
    cal=pm.calibrate()
    DA_dir=pm.DA_dir()
    os.system("source "+pm.DA_dir()+"/src/oneday_sim.sh "+yyyy+" "+mm+" "+dd+" "+ens_num+" "+dir2
    +" "+looptype+" "+str(cpunums)+" "+str(run_name)+" "+str(exp_dir)+" "+str(mapname)+" "+str(cal)
    +" "+DA_dir)

    return 0
########################### # modified to run paralle @Menaka 
def copy_corrupted_sfcelv(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1] 
    dd   = inputlist[2]
    num  = inputlist[3]
    numch='%03d'%num
    fname="./CaMa_out/"+yyyy+mm+dd+"C"+numch+"/sfcelv"+yyyy+".bin"
    os.system("cp "+fname+" ./assim_out/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin")
    return 0
########################### # modified not calculate restart again/ no chage in WSE in corrupted @Menaka
def copy_corrupted_restart(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    num  = inputlist[3]
    nxt_day = datetime.date(int(yyyy),int(mm),int(dd)) + datetime.timedelta(days=1)
    n_yyyy='%04d' % (nxt_day.year)
    n_mm='%02d' % (nxt_day.month)
    n_dd='%02d' % (nxt_day.day)
    numch='%03d'%num
    fname="./CaMa_out/"+yyyy+mm+dd+"C"+numch+"/restart"+n_yyyy+n_mm+n_dd+".bin"
    #os.system("cp "+fname+" ./CaMa_in/restart/open/restart"+n_yyyy+n_mm+n_dd+"C"+numch+".bin")
    copy_stoonly(fname,"./CaMa_in/restart/open/restart"+n_yyyy+n_mm+n_dd+"C"+numch+".bin")
    print ("copy restart",n_yyyy,n_mm,n_dd,"C"+numch)
    return 0
###########################
def copy_stoonly(iname,oname): # for CaMa_Flood v395b
    org=np.fromfile(iname,np.float32).reshape(6,-1)
    org[0:2].tofile(oname)
    return 0
###########################
def data_assim(yyyy,mm,dd): # new data assimilation function (2020/05/18)
    #print '%02d'%(nxt_day.day)
    parallels="%d"%(pm.para_nums()*pm.cpu_nums())
    os.environ['OMP_NUM_THREADS']=parallels
    #os.system("export $OMP_NUM_THREADS=%d"%(pm.para_nums()*pm.cpu_nums()))
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    print (pm.ens_mem(pm.mode()))
    thisday=datetime.date(int(yyyy),int(mm),int(dd))
    nxt_day=thisday+datetime.timedelta(days=1)
    nyear=nxt_day.year
    nmon=nxt_day.month
    nday=nxt_day.day
    nxtyyyymmdd="%04d%02d%02d"%(nyear,nmon,nday)
    os.system(pm.DA_dir()+"/src/data_assim "+yyyy+mm+dd+" "+pm.mapname()+" "\
    +str(pm.patch_size())+" "+str(pm.ens_mem(pm.mode()))+" "+nxtyyyymmdd+" "+pm.CaMa_dir()\
    +" "+str(pm.thersold())+" "+exp_dir+" "+pm.DA_dir()+" "+pm.patch_dir()+" "\
    +str(pm.patch_name())+" "+pm.HydroWeb_dir()+" "+str(pm.rho())+" "+str(pm.sigma_b())\
    +" "+str(pm.conflag())+" "+(pm.calibrate()))
    return 0
###########################
def make_init_storge():
    bef_yyyy='%04d' % (pm.spinup_end_year())
    bef_mm='%02d' % (pm.spinup_end_month())
    bef_dd='%02d' % (pm.spinup_end_date())
    os.system("./src/make_nonassim_init ./CaMa_out/spinup_open/storge"+str(pm.spinup_end_year())+".bin "+"./assim_out/nonassim/open/nonasmC"+bef_yyyy+bef_mm+bef_dd+".bin")
    os.system("./src/make_nonassim_init ./CaMa_out/spinup_open/storge"+str(pm.spinup_end_year())+".bin "+"./assim_out/nonassim/assim/nonasmA"+bef_yyyy+bef_mm+bef_dd+".bin")
    return 0
###########################
def make_initial_restart(): # updated the name
    start_year,start_month,start_date=pm.starttime()
    yyyy="%04d"%(start_year)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    #spinup_true="%04d%2d%02dT000"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date()) 
    #os.system("cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    #copy_stoonly(exp_dir+"/CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin",exp_dir+"/CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")

    #print "cp "+exp_dir+"/CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin  "+exp_dir+"/CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin"
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        spinup_open="%04d%2d%02dC%03d"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date(),num) 
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
        copy_stoonly(exp_dir+"CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin",exp_dir+"CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        copy_stoonly(exp_dir+"CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin",exp_dir+"CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
###########################
def make_initial_restart_one(): # updated the name
    # copy restartyyyymmddC001 as restart for all simulations 
    start_year,start_month,start_date=pm.starttime()
    yyyy="%04d"%(start_year)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    spinup_true="%04d%2d%02dT000"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date())
    #os.system("cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    copy_stoonly("./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    print "cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin"
    spinup_open=spinup_true
    #spinup_open="%04d%2d%02dC%03d"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date(),1)
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
###########################