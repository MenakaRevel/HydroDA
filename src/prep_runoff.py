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
# #########################################
#
# Prepare input runoff fields fo CaMa-Flood
#
# #########################################
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
def copy_runoff(inputlist): #do it parallel
    iname=inputlist[0]
    oname=inputlist[1]
    distopen=float(inputlist[2])
    runoff=np.fromfile(iname,np.float32)
    runoff=(np.ma.masked_equal(runoff,-9999.0)*distopen).filled(-9999.0)
    runoff.tofile(oname)
    return 0
###########################
def cal_monthly_mean_ens(ens_num): #start_year,end_year,ens_num,months=24):
    # calc monthly mean value for two years
    start_year=pm.spinup_end_year()
    end_year,end_month,end_date=pm.starttime()
    #start_year=pm.start_year() #inputlist[0]
    #end_year=pm.end_year() #inputlist[1]
    #ens_num=inputlist[2]
    months=24 #inputlist[3]
    #threshold=inputlist[4]
    runname=pm.runname(pm.mode())
    if runname=="E2O":
        nx=1440
        ny=720
        threshold=0.1
    else:
        nx=360
        ny=180
        threshold=1.0e-8
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    #threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    ens_char="%03d"%(ens_num)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile(pm.DA_dir()+"/inp/"+runname+"/Roff/Roff__"+str(yyyy)+str(mm)+str(dd)+ens_char+".one",np.float32).reshape([ny,nx])
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile("./CaMa_in/"+runname+"/mean_month/mean_"+ychar+mchar+ens_char+".bin")
###########################
def cal_monthly_total(start_year,end_year,months=24,threshold=0.1):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    true_run=pm.true_run(pm.mode())
    if runname=="ELSE_KIM2009":
        nx,ny=360,180
        prefix="Roff____"
        suffix=".one"
    if runname=="E2O":
        nx,ny=1440,720
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    if runname=="ERA20CM":
        nx,ny=360,180
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    mkdir("./CaMa_in/"+runname+"/total_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    #os.system("rm -Rf ./CaMa_in/"+runname+"/total_month/*")
    #threshold=0.1
    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        #count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/"+prefix+str(yyyy)+str(mm)+str(dd)+suffix,np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            #count=count+(roff>threshold)
        roff_total=roff_mon #/(count+1e-20)
        roff_total=roff_total.astype(np.float32)
        roff_total=roff_total+threshold
        roff_total.tofile("./CaMa_in/"+runname+"/total_month/total_"+ychar+mchar+".bin")
###########################
def cal_monthly_mean(start_year,end_year,months=24):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    true_run=pm.true_run(pm.mode())
    if runname=="ELSE_KIM2009":
        nx,ny=360,180
        prefix="Roff____"
        suffix=".one"
    if runname=="E2O":
        nx,ny=1440,720
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    if runname=="ERA20CM":
        nx,ny=360,180
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/"+prefix+str(yyyy)+str(mm)+str(dd)+suffix,np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile("./CaMa_in/"+runname+"/mean_month/mean_"+ychar+mchar+".bin")
###########################
def prepare_input_old():
    # spinup start_year
    # simulation end_year
    start_year=pm.spinup_end_year()
    end_year,end_month,end_date=pm.endtime()
    #start_year=pm.start_year()
    #end_year=pm.end_year()
    start_dt=datetime.date(start_year,1,1)
    last_dt=datetime.date(end_year,end_month,end_date)
    start=0
    last=int((last_dt-start_dt).days)
    #--------------
    # E2O
    if pm.input()=="E2O": # Earth2Observe
        #print "E2O"
        distopen=pm.distopen(1)
        diststd=pm.diststd(1)
        true_run=pm.true_run(1) # for true ensemble
        runname=pm.runname(1) # get runoff name
        # make courrpted/assimilated runoff
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            #print "Roff_CORR available"
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR")
        #print "E2O/Roff"
        # # # dist std for ensembles
        # # distopen_ranges={}
        # # #open_list=np.setdiff1d(np.arange(1,7+1),[true_run])
        # # open_list=np.arange(1,7+1)
        # # random_runs=random.sample(open_list,k=2)
        # # mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/runoff_error")
        # # f=open(pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/runoff_error/ensemble.txt","w")
        # # for runens in open_list:
        # #     ens_char="%03d"%(runens)
        # #     diststd_num=3
        # #     distopen_range=rd.normal(1,diststd,diststd_num)
        # #     distopen_range=np.sort(distopen_range)
        # #     distopen_range=distopen_range.astype(np.float32)
        # #     distopen_ranges[ens_char]=distopen_range#[0.25,1.00,1.25]
        # #     line="%s  %3.2f  %3.2f  %3.2f\n"%(ens_char,distopen_range[0],distopen_range[1],distopen_range[2])
        # #     f.write(line)
        # # f.close()
        # # #print "L1141"
        # # inputlist=[]
        # # for day in np.arange(start,last):
        # #     target_dt=start_dt+datetime.timedelta(days=day)
        # #     yyyy='%04d' % (target_dt.year)
        # #     mm='%02d' % (target_dt.month)
        # #     dd='%02d' % (target_dt.day)
        # #     ens_num=1
        # #     for runens in open_list: #np.arange(1,7+1):
        # #         #print runens
        # #         #if runens!=true_run:
        # #         run_num="%03d"%(runens)
        # #         iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff__"+yyyy+mm+dd+run_num+".one"
        # #         distopens=distopen_ranges[run_num]
        # #         for dist in distopens:
        # #             ens_char="C%03d"%(ens_num)
        # #             oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
        # #             inputlist.append([iname,oname,str(abs(dist))])
        # #             ens_num=ens_num+1
        distopen_ranges={}
        # fname="random_ensemble_E2O.txt"
        fname="random_ensemble_E2O_49.txt"
        with open(pm.DA_dir()+"/dat/"+fname,"r") as f:
            lines=f.readlines()
        for line in lines:
            line   = filter(None,re.split(" ", line))
            # print line
            runens = int(line[0])
            rndnum = float(line[1].strip("\n"))
            ens_char="C%03d"%(runens)
            distopen_ranges[ens_char]=rndnum
        #print "L1141"
        inputlist=[]
        open_list=np.arange(1,7+1)
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_num=1
            for runens in open_list: #np.arange(1,7+1):
                #print runens
                #if runens!=true_run:
                run_num="%03d"%(runens)
                iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff__"+yyyy+mm+dd+run_num+".one"
                # distopens=distopen_ranges[run_num]
                # for dist in distopens:
                for _ in np.arange(1,7+1):
                    ens_char="C%03d"%(ens_num)
                    dist=distopen_ranges[ens_char]
                    oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                    inputlist.append([iname,oname,str(abs(dist))])
                    ens_num=ens_num+1
        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()
    #---------
    # ERA20CM
    if pm.input()=="ERA20CM": #ECMWF ERA20CM
        distopen=pm.distopen(2)
        diststd=pm.diststd(2)
        true_run=pm.true_run(2) # for true ensemble
        runname=pm.runname(2) # get runoff name
        # make courrpted runoff
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR")
        #--
        # dist std for ensembles
        distopen_ranges={}
        #open_list=np.setdiff1d(np.arange(1,10+1),[true_run])
        open_list=np.arange(1,10+1)
        random_runs=random.sample(open_list,k=2)
        for runens in open_list:
            ens_char="%03d"%(runens)
            diststd_num=2
            distopen_range=rd.normal(1,diststd,diststd_num)
            distopen_range=np.sort(distopen_range)
            distopen_range=distopen_range.astype(np.float32)
            distopen_ranges[ens_char]=distopen_range#[0.25,1.00,1.25]
        #
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_num=1
            for runens in open_list: #np.arange(1,10+1):
                run_num="%03d"%(runens)
                iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff__"+yyyy+mm+dd+run_num+".one"
                distopens=distopen_ranges[run_num]
                for dist in distopens:
                    ens_char="C%03d"%(ens_num)
                    oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                    inputlist.append([iname,oname,str(dist)])
                    ens_num=ens_num+1
        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()
    #---------
    # VIC BC
    if pm.input()=="VIC_BC": #VIC BC
        nXX,nYY=1440,720
        distopen=1.0 #pm.distopen(3)
        diststd=0.25 #pm.diststd(3)
        true_run=pm.true_run(3) # for true ensemble
        runname=pm.runname(3) # get runoff name
        # make courrpted runoff
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR")
        #--
        # dist std for ensembles
        #distopen_ranges={}
        # for 12 months
        distopen_range=np.zeros([12,pm.ens_mem(),nYY,nXX],np.float32)
        # fname="../../dat/std_runoff_E2O_1980-2000.bin"
        # #print fname
        # std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
        # fname="../../dat/mean_runoff_E2O_1980-2000.bin"
        # #print "L1481",fname
        # mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
        # #print mean_runoff
        # std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
        # mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
        for mon in range(1,12+1): # for 12 months
            mm="%02d"%(mon)
            fname=pm.DA_dir()+"/dat/std_month_runoff_E2O_"+mm+".bin"
            #print fname
            std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            fname=pm.DA_dir()+"/dat/mean_month_runoff_E2O_"+mm+".bin"
            #print "L1481",fname
            mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            #print mean_runoff
            std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
            mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
            for iXX in range(nXX):
                for iYY in range(nYY):
                    #distopen_range[:,iYY,iXX]=np.sort(rd.normal(distopen,std_runoff[iYY,iXX],pm.ens_mem()))
                    #Log-normal model
                    #sk=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
                    sk=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
                    beta=0.0
                    #E=std_runoff[iYY,iXX]/(mean_runoff[iYY,iXX]+1.0e-20)
                    E=diststd
                    #distopen_range[mon,:,iYY,iXX]=((1+beta)/math.sqrt(E**2+1))*np.exp(math.sqrt(math.log(E**2+1))*sk)
                    distopen_range[mon-1,:,iYY,iXX]=np.sort(rd.normal(distopen,E,pm.ens_mem()))
                    # distopen_range[mon-1,:,iYY,iXX]=rd.normal(distopen,E,pm.ens_mem())
                    #distopen_range[:,iYY,iXX]=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
        #----------
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            mon=int(target_dt.month)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one"
            #print iname
            roff=np.fromfile(iname,np.float32).reshape(nYY,nXX)
            #roff=np.ones([nYY,nXX],np.float32)*-9999.0
            fname=pm.DA_dir()+"/dat/std_runoff_E2O_1980-2000.bin"
            #print fname
            std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            fname=pm.DA_dir()+"/dat/mean_runoff_E2O_1980-2000.bin"
            #print "L1481",fname
            mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            #print mean_runoff
            std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
            mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
            for ens_num in np.arange(pm.ens_mem()):
                ens_char="C%03d"%(ens_num+1)
                roffc=roff+distopen_range[mon-1,ens_num,:,:]*mean_runoff #*std_runoff #
                roffc=roff*distopen_range[mon-1,ens_num,:,:] #*mean_runoff
                oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                roffc.tofile(oname)
    #--------------
    # -25% biased runoff experiment
    if pm.input()=="ELSE": # Kim 2009/E2O/ERA20CM
        #print "-25% biased runoff experiment", pm.true_run(3), pm.runname(3)
        distopen=pm.distopen(3)
        diststd=pm.diststd(3)
        true_run=pm.true_run(3) # for true ensemble
        runname=pm.runname(3) # get runoff name
        # copy for TRUE simulation
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE")
        #print true_run
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            true_char="%03d"%(true_run)
            if runname=="ELSE_KIM2009":
                nx,ny=360,180
                prefix="Roff____"
                suffix=".one"
            if runname=="E2O":
                nx,ny=1440,720
                prefix="Roff__"
                suffix="%03d.one"%(true_run)
            if runname=="ERA20CM":
                nx,ny=360,180
                prefix="Roff__"
                suffix="%03d.one"%(true_run)
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/"+prefix+yyyy+mm+dd+suffix
            oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate mothly mean
        cal_monthly_mean(start_year,end_year,24)

        # calculate monthly total
        cal_monthly_total(start_year,end_year,24)

        # make courrpted runoff
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR")

        inputlist=[]
        std=rd.normal(0,diststd,pm.ens_mem())
        std=std.astype(np.float32)
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # make random values
            #std=np.fromfile("./CaMa_out/randlist.bin",np.float32)
            #std=rd.normal(0,diststd,pm.ens_mem())
            #std=std.astype(np.float32)
            #print std
            ifile=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+"T000.one",np.float32).reshape(ny,nx)
            roff_mean=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            roff_total=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            for ens in np.arange(1,pm.ens_mem()+1):
                ens_char="C%03d"%(ens)
                #print pm.distopen(),std[ens-1]
                ofile=ifile*distopen + roff_mean*std[ens-1]#*10.0
                #ofile=ifile*(distopen + std[ens-1])
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile.astype(np.float32)
                ofile.tofile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    #--------------
    # blind runoff experiment
    if pm.input()=="ELSE_KIM2009": # differnt yeaer
        #distopen=1.0 #0.75 #pm.distopen()
        #diststd=0.5  #pm.diststd()
        distopen=pm.distopen(4)
        diststd=pm.diststd(4)
        true_run=pm.true_run(4) # for true ensemble
        runname=pm.runname(4) # get runoff name  
        # copy for TRUE simulation
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one"
            oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate mothly mean
        cal_monthly_mean(start_year,end_year,24)

        # calculate monthly total
        cal_monthly_total(start_year,end_year,24)

        # make courrpted runoff
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # make random values
            #std=np.fromfile("./CaMa_out/randlist.bin",np.float32)
            std=rd.normal(0,diststd,pm.ens_mem())
            std=std.astype(np.float32)
            #print std
            ifile=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one",np.float32).reshape(ny,nx)
            roff_mean=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            roff_total=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            for ens in np.arange(1,pm.ens_mem()+1):
                ens_char="C%03d"%(ens)
                #print pm.distopen(),std[ens-1]
                ofile=ifile*distopen + roff_mean*std[ens-1]*10.0
                #ofile=ifile*(distopen + std[ens-1])
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile.astype(np.float32)
                ofile.tofile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    #---------
    # ERA5
    if pm.input()=="ERA5": #ERA5-Land
        nXX,nYY=3600,1800
        distopen=pm.distopen(5)
        diststd=pm.diststd(5)
        true_run=pm.true_run(3) # for true ensemble
        runname=pm.runname(3) # get runoff name
        # make courrpted runoff
        if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR")
        #--
        # dist std for ensembles
        #distopen_ranges={}
        # for 12 months
        distopen_range=np.zeros([12,pm.ens_mem(),nYY,nXX],np.float32)
        # fname="../../dat/std_runoff_E2O_1980-2000.bin"
        # #print fname
        # std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
        # fname="../../dat/mean_runoff_E2O_1980-2000.bin"
        # #print "L1481",fname
        # mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
        # #print mean_runoff
        # std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
        # mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
        for mon in range(1,12+1): # for 12 months
            mm="%02d"%(mon)
            fname=pm.DA_dir()+"/dat/std_month_runoff_E2O_"+mm+".bin"
            #print fname
            std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            fname=pm.DA_dir()+"/dat/mean_month_runoff_E2O_"+mm+".bin"
            #print "L1481",fname
            mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            #print mean_runoff
            std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
            mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
            for iXX in range(nXX):
                for iYY in range(nYY):
                    #distopen_range[:,iYY,iXX]=np.sort(rd.normal(distopen,std_runoff[iYY,iXX],pm.ens_mem()))
                    #Log-normal model
                    #sk=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
                    sk=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
                    beta=0.0
                    #E=std_runoff[iYY,iXX]/(mean_runoff[iYY,iXX]+1.0e-20)
                    E=diststd
                    #distopen_range[mon,:,iYY,iXX]=((1+beta)/math.sqrt(E**2+1))*np.exp(math.sqrt(math.log(E**2+1))*sk)
                    distopen_range[mon-1,:,iYY,iXX]=np.sort(rd.normal(distopen,E,pm.ens_mem()))
                    # distopen_range[mon-1,:,iYY,iXX]=rd.normal(distopen,E,pm.ens_mem())
                    #distopen_range[:,iYY,iXX]=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
        #----------
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            mon=int(target_dt.month)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".sixmin"
            #print iname
            roff=np.fromfile(iname,np.float32).reshape(nYY,nXX)
            #roff=np.ones([nYY,nXX],np.float32)*-9999.0
            fname=pm.DA_dir()+"/dat/std_runoff_E2O_1980-2000.bin"
            #print fname
            std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            fname=pm.DA_dir()+"/dat/mean_runoff_E2O_1980-2000.bin"
            #print "L1481",fname
            mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            #print mean_runoff
            std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
            mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
            for ens_num in np.arange(pm.ens_mem()):
                ens_char="C%03d"%(ens_num+1)
                roffc=roff+distopen_range[mon-1,ens_num,:,:]*mean_runoff #*std_runoff #
                roffc=roff*distopen_range[mon-1,ens_num,:,:] #*mean_runoff
                oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                roffc.tofile(oname)
    return 0
######################################
def prepare_input():
    """
    link to runoff data set
    """
    if os.path.islink("./CaMa_in/"+pm.runname(pm.mode())):
        os.system("rm -r ./CaMa_in/"+pm.runname(pm.mode()))
    os.system("ln -sf "+pm.runoff_dir()+" ./CaMa_in/")
    return 0
######################################
if __name__ == "__main__":
    print ("prepare runoff input")
    prepare_input()