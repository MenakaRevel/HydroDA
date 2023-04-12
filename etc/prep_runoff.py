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
import scipy.signal

# #external python codes
# sys.path.append("../gosh")
# import params as pm
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
def cal_monthly_total(runname,start_year,end_year,indir,outdir,months=24,threshold=0.1):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    true_run=pm.true_run(pm.mode())
    nx,ny,prefix,suffix=get_runoff_metadata(runname)
    # if runname=="ELSE_KIM2009":
    #     nx,ny=360,180
    #     prefix="Roff____"
    #     suffix=".one"
    # if runname=="E2O":
    #     nx,ny=1440,720
    #     prefix="Roff__"
    #     suffix="%03d.one"%(true_run)
    # if runname=="ERA20CM":
    #     nx,ny=360,180
    #     prefix="Roff__"
    #     suffix="%03d.one"%(true_run)
    
    mkdir(outdir+"/total_month")
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
            roff=np.fromfile(indir+"/"+prefix+str(yyyy)+str(mm)+str(dd)+suffix,np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            #count=count+(roff>threshold)
        roff_total=roff_mon #/(count+1e-20)
        roff_total=roff_total.astype(np.float32)
        roff_total=roff_total+threshold
        roff_total.tofile(outdir+"/"+runname0+"/total_month/total_"+ychar+mchar+".bin")
###########################
def cal_monthly_mean(runname0,syear,eyear,indir,outdir,months=24):
    # calc monthly mean value for two years
    # runname0=pm.runname(pm.mode())
    true_run=1 #pm.true_run(pm.mode())
    nx,ny,prefix,suffix=get_runoff_metadata(runname0)
    # if runname0=="ELSE_KIM2009":
    #     nx,ny=360,180
    #     prefix="Roff____"
    #     suffix=".one"
    #     runname=runname0
    # if runname0=="E2O":
    #     nx,ny=1440,720
    #     prefix="Roff__"
    #     suffix="%03d.one"%(true_run)
    #     runname=runname0
    # if runname0=="ERA20CM":
    #     nx,ny=360,180
    #     prefix="Roff__"
    #     suffix="%03d.one"%(true_run)
    #     runname=runname0
    # if runname0=="ECMWF":
    #     nx,ny=1440,720
    #     prefix="Roff__"
    #     suffix="003.one"
    #     runname="E2O"
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir(outdir+"/"+runname0+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    threshold=0.1
    #--
    start_dt=datetime.date(syear,1,1)
    end_dt=datetime.date(eyear,12,31)
    for month in np.arange(months):
        ynow=int(syear+int(month/12))
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
            roff=np.fromfile(indir+"/"+prefix+str(yyyy)+str(mm)+str(dd)+suffix,np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile(outdir+"/"+runname0+"/mean_month/mean_"+ychar+mchar+".bin")
        # print (outdir+"/"+runname0+"/mean_month/mean_"+ychar+mchar+".bin")
    return 0
###########################
def prepare_input():
    # spinup start_year
    # simulation end_year
    # ISMIP3a
    start_year=pm.start_year()
    end_year=pm.end_year()
    #start_year=pm.start_year()
    #end_year=pm.end_year()
    start_dt=datetime.date(start_year,1,1)
    last_dt=datetime.date(end_year,12,31)
    start=0
    last=int((last_dt-start_dt).days)
    #--------------
    # E2O
    if pm.input()=="E2O": # Earth2Observe
        #print "E2O"
        # distopen=pm.distopen(1)
        # diststd=pm.diststd(1)
        # true_run=pm.true_run(1) # for true ensemble
        runname=pm.runname(1) # get runoff name
        # make courrpted/assimilated runoff
        if(len("./CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0:
            #print "Roff_CORR available"
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff")
        #print "E2O/Roff"
        # dist std for ensembles
        distopen_ranges={}
        #open_list=np.setdiff1d(np.arange(1,7+1),[true_run])
        # open_list=np.arange(1,7+1)
        # random_runs=random.sample(open_list,k=2)
        # mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/runoff_error")
        fname="./ensemble_E2O_49.txt"
        with open(fname,"r") as f:
            lines=f.readlines()
        for line in lines:
            line   = filter(None,re.split(" ", line))
            # print line
            runens = int(line[0])
            rndnum = float(line[1].strip("\n"))
            ens_char="%03d"%(runens)
            distopen_ranges[ens_char]=rndnum
        # for runens in open_list:
        #     ens_char="%03d"%(runens)
        #     diststd_num=3
        #     distopen_range=rd.normal(1,diststd,diststd_num)
        #     distopen_range=np.sort(distopen_range)
        #     distopen_range=distopen_range.astype(np.float32)
        #     distopen_ranges[ens_char]=distopen_range#[0.25,1.00,1.25]
        #     line="%s  %3.2f  %3.2f  %3.2f\n"%(ens_char,distopen_range[0],distopen_range[1],distopen_range[2])
        #     f.write(line)
        # f.close()
        #print "L1141"
        inputlist=[]
        for day in np.arange(start,last+1):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_num=1
            for run_num in np.arange(1,7+1): #np.arange(1,7+1):
                #print runens
                #if runens!=true_run:
                run_char="%03d"%(run_num)
                iname=pm.inputdir()+"/"+runname+"/Roff/Roff__"+yyyy+mm+dd+run_char+".one"
                for _ in np.arange(1,7+1):
                    ens_char="%03d"%(ens_num)
                    dist=distopen_ranges[ens_char]
                    # for dist in distopens:
                    # ens_char="C%03d"%(ens_num)
                    # print run_char, ens_char
                    oname="./CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
                    print iname, oname, str(abs(dist))
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
    #---------
    # ECMWF
    if pm.input()=="ECMWF": #ECMWF
        print pm.input()
        nXX,nYY=1440,720
        distopen=pm.distopen()
        diststd=pm.diststd()
        # true_run=pm.true_run(3) # for true ensemble
        runname=pm.runname(5) # get runoff name
        # make courrpted runoff
        if(len(glob.glob("./CaMa_in/"+runname+"/Roff/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff")

        # calculate mothly mean
        cal_monthly_mean(start_year,end_year,(end_year-start_year)*12)

        # calculate monthly total
        # cal_monthly_total(start_year,end_year,24)

        inputlist=[]
        std=np.random.normal(0,diststd,pm.ens_mem())
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
            # iname=pm.inputdir()+"/"+runname+"/Roff/Roff__"+yyyy+mm+dd+run_char+".one"
            ifile=np.fromfile(pm.inputdir()+"/E2O/Roff/Roff__"+yyyy+mm+dd+"003.one",np.float32).reshape(nYY,nXX)
            roff_mean=np.fromfile("./CaMa_in/"+runname+"/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(nYY,nXX)
            # roff_total=np.fromfile(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            for ens in np.arange(1,pm.ens_mem()+1):
                ens_char="%03d"%(ens)
                print (pm.distopen(),std[ens-1])
                ofile=ifile*distopen + roff_mean*std[ens-1]#*10.0
                #ofile=ifile*(distopen + std[ens-1])
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile.astype(np.float32)
                ofile.tofile("./CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one")
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
    #--------------
    # ERA5
    if pm.input()=="ERA5": # ERA5
        global rand
        nx,ny=3600,1800
        #print "E2O"
        # distopen=pm.distopen(1)
        # diststd=pm.diststd(1)
        # true_run=pm.true_run(1) # for true ensemble
        runname=pm.runname(1) # get runoff name
        # make courrpted/assimilated runoff
        if(len("./CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0:
            #print "Roff_CORR available"
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff")
        #print "E2O/Roff"
        # dist std for ensembles
        distopen_ranges={}
        nx_char="%d"%(nx)
        ny_char="%d"%(ny)
        #open_list=np.setdiff1d(np.arange(1,7+1),[true_run])
        # open_list=np.arange(1,7+1)
        # random_runs=random.sample(open_list,k=2)
        # mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/runoff_error")
        rivnum=np.fromfile("./dat/rivnum_glb_06min.bin",np.float32).reshape(ny,nx)
        rand=np.ones([ny,nx,pm.ens_mem()])
        rd=np.random.normal(1.0,0.0025,pm.ens_mem())
        for ens in np.arange(0,pm.ens_mem()):
            rand[:,:,ens]=rd[ens]
        with open("./random_subbasin_id_E2O_std.txt","r") as f:
            lines=f.readlines()
        for line in lines:
            line   = filter(None,re.split(" ", line))
            # print line
            id     = float(line[0])
            std    = float(line[1])
            runens = int(line[1])
            indice = np.where(rivnum==id)
            for ens in np.arange(1,pm.ens_mem()+1):
                randnum = float(line[ens+1].strip("\n"))
                rand[indice[0],indice[1],ens-1] = randnum
            # ens_char="%03d"%(runens)
            # distopen_ranges[ens_char]=rndnum
        # for runens in open_list:
        #     ens_char="%03d"%(runens)
        #     diststd_num=3
        #     distopen_range=rd.normal(1,diststd,diststd_num)
        #     distopen_range=np.sort(distopen_range)
        #     distopen_range=distopen_range.astype(np.float32)
        #     distopen_ranges[ens_char]=distopen_range#[0.25,1.00,1.25]
        #     line="%s  %3.2f  %3.2f  %3.2f\n"%(ens_char,distopen_range[0],distopen_range[1],distopen_range[2])
        #     f.write(line)
        # f.close()
        #print "L1141"
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # ens_num=1
            for ens in np.arange(1,pm.ens_mem()+1): #np.arange(1,7+1):
                iname=pm.inputdir()+runname+"/bin/Roff____"+yyyy+mm+dd+".sixmin"
                ens_char="%03d"%(ens)
                oname="./CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
                inputlist.append([iname,oname,ens_char,nx_char,ny_char])
                # ens_num=ens_num+1
        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff_subbasin,inputlist)
        p.terminate()
    #---------
    return 0
######################################
def prep_runoff_ensemble(distopen,diststd,ne,runname,rundir,
                        outdir,syear=2001,eyear=2010,method="simple",
                        beta=0.0,E=0.3,alpha=0.993,tau_s=50.0):
    """
    runname = name of the runoff [e.g., VIC BC, E2O]
    rundir  = runoff directory
    outdir  = directory to save output
    nx, ny  = spatial resolution
    #=========================================================================
    A function to create runoff pertubations
    Several methods propsed
    1. simple - random normal distribution spatial constant
    2. spatial - spatially correlated random field
    3. lognorm - spatially correlated random fields by log normal distribution as proposed by Nijssen and Lettenmaier (2004)
    """
    #=================================
    # get runoff metadata
    nx,ny,prefix,sufix=get_runoff_metadata(runname)
    # make directories
    mkdir(outdir+"/"+runname)
    mkdir(outdir+"/"+runname+"/Roff")
    # create necessary varibales
    start_dt=datetime.date(syear,1,1)
    end_dt=datetime.date(eyear,12,31)
    start=0
    last=(end_dt-start_dt).days + 1
    #=================================
    if method=="simple":
        #distopen=1.0 #0.75 #pm.distopen()
        #diststd=0.5  #pm.diststd()
        # distopen=pm.distopen(4)
        # diststd=pm.diststd(4)
        # true_run=pm.true_run(4) # for true ensemble
        # runname=pm.runname(4) # get runoff name
        # ne=20 #pm.ens_mem() # nuumber of pertubations
        # # get runoff metadata
        # nx,ny,prefix,sufix=get_runoff_metadata(runname)
        # # make directories
        # mkdir(outdir+"/"+runname)
        # mkdir(outdir+"/"+runname+"/Roff")
        # get random values

        # # copy for TRUE simulation
        # if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff*"))!=0):
        #     pass
        # # true input file is not ready
        # need preparation
        # make directories
        # mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE")
        # start_dt=datetime.date(syear,1,1)
        # end_dt=datetime.date(eyear,12,31)
        # start=0
        # last=(end_dt-start_dt).days
        # inputlist=[]
        # for day in np.arange(start,last):
        #     target_dt=start_dt+datetime.timedelta(days=day)
        #     yyyy='%04d' % (target_dt.year)
        #     mm='%02d' % (target_dt.month)
        #     dd='%02d' % (target_dt.day)
        #     ens_char="T000"
        #     iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one"
        #     oname=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
        #     inputlist.append([iname,oname,"1.00"])

        # # do parallel
        # p=Pool(pm.para_nums())
        # p.map(copy_runoff,inputlist)
        # p.terminate()

        # calculate mothly mean
        cal_monthly_mean(runname,syear,eyear,rundir,outdir,12*(eyear-syear+1))

        # calculate monthly total
        # cal_monthly_total(runname,syear,eyear,rundir,outdir,24)

        # make runoff pertubation
        # if(len(glob.glob(pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
        #     pass
        # corrupted input file is not ready
        # need preparation
        
        #====================
        # # get runoff metadata
        # nx,ny,prefix,suffix=get_runoff_metadata(runname)
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # make random values
            #std=np.fromfile("./CaMa_out/randlist.bin",np.float32)
            std=rd.normal(0,diststd,ne)
            std=std.astype(np.float32)
            #print std
            # names
            ifile=rundir+"/"+prefix+yyyy+mm+dd+suffix
            roff_mean=outdir+"/"+runname+"/mean_month/mean_"+yyyy+mm+".bin"
            # ifile=np.fromfile(rundir+"/"+prefix+yyyy+mm+dd+suffix,np.float32).reshape(ny,nx)
            # roff_mean=np.fromfile(outdir+"/"+runname+"/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            # roff_total=np.fromfile(outdir+"/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            #===========================
            for ens in np.arange(1,ne+1):
                ens_char="%03d"%(ens)
                #print pm.distopen(),std[ens-1]
                # ofile=ifile*distopen + roff_mean*std[ens-1]
                #ofile=ifile*(distopen + std[ens-1])
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile=outdir+"/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
                # ofile.astype(np.float32)
                # ofile.tofile(outdir+"/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")
                inputlist.append([yyyy,mm,dd,ens_char,str(distopen),str(std[ens-1]),ifile,roff_mean,str(nx),str(ny),ofile])
        #==================================
        # do parallel
        para=20 #pm.para_nums()
        p=Pool(para)
        p.map(make_runoff,inputlist)
        p.terminate()
        # map(make_runoff,inputlist)
    #=================================
    if method=="lognormal":
        # spatially distibuted lognormal
        # get runoff metadata
        nx,ny,prefix,sufix=get_runoff_metadata(runname)
        # make directories
        mkdir(outdir+"/"+runname)
        mkdir(outdir+"/"+runname+"/Roff")
        #====================
        # get runoff metadata
        nx,ny,prefix,suffix=get_runoff_metadata(runname)
        #====================
        # alpha=1.0-(1.0/(tau_t+1e-20))
        # beta=val_beta()
        # E=val_E()
        #===========================
        s=spatially_correlated_random(nx,ny)
        rand=np.sort(np.abs(np.random.normal(0,1,[ne])))
        inputlist=[]
        for ens in np.arange(1,ne+1):
            ens_char="%03d"%(ens)
            yyyy1="%04d"%(syear)
            yyyy2="%04d"%(eyear)
            # print (yyyy1,yyyy2,ens_char,runname,str(alpha),str(beta),str(E),rundir,outdir)
            inputlist.append([yyyy1,yyyy2,ens_char,runname,str(alpha),str(beta),str(E),str(tau_s),rundir,outdir])
        #==================================
        # do parallel
        para=20 #pm.para_nums()
        p=Pool(para)
        p.map(make_runoff_temporal,inputlist)
        p.terminate()

            # for day in np.arange(start,last):
            #     target_dt=start_dt+datetime.timedelta(days=day)
            #     yyyy='%04d' % (target_dt.year)
            #     mm='%02d' % (target_dt.month)
            #     dd='%02d' % (target_dt.day)
            #     # get spatially correlated normal distributiion
            #     w=spatially_correlated_random(nx,ny)
            #     s=alpha*s+math.sqrt(1-alpha**2)*w
            #     Rcoff=((1+beta)/math.sqrt(E**2+1))*math.exp(math.sqrt(E**2+1))*s
            #     # ifile=rundir+"/"+prefix+yyyy+mm+dd+suffix
            #     ifile=np.fromfile(rundir+"/"+prefix+yyyy+mm+dd+suffix,np.float32).reshape(ny,nx)
            #     ofile=outdir+"/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
            #     Rc=Rcoff*ifile
            #     Rc.astype(np.float32)
            #     Rc.tofile(ofile)
    #=================================
    if method=="normal":
        # spatially distibuted lognormal
        # get runoff metadata
        nx,ny,prefix,sufix=get_runoff_metadata(runname)
        # make directories
        mkdir(outdir+"/"+runname)
        mkdir(outdir+"/"+runname+"/Roff")
        #====================
        # get runoff metadata
        nx,ny,prefix,suffix=get_runoff_metadata(runname)
        #====================
        #===========================
        # s=spatially_correlated_random_init(nx,ny)
        # rand=np.sort(np.abs(np.random.normal(0,1,[ne])))
        inputlist=[]
        for ens in np.arange(1,ne+1):
            ens_char="%03d"%(ens)
            yyyy1="%04d"%(syear)
            yyyy2="%04d"%(eyear)
            # print (yyyy1,yyyy2,ens_char,runname,str(alpha),str(beta),str(E),rundir,outdir)
            inputlist.append([yyyy1,yyyy2,ens_char,runname,str(alpha),str(beta),str(E),str(tau_s),rundir,outdir])
        #==================================
        # do parallel
        para=20 #pm.para_nums()
        p=Pool(para)
        p.map(make_runoff_normal,inputlist)
        p.terminate()
    return 0
######################################
def spatially_correlated_random(nx,ny,tau=50.0):
    """
    Spatially correlated random variable 
    sptial demensions are [nx,ny]
    tau_s = decorrelation length
    """
    correlation_scale = tau
    x = np.arange(-correlation_scale, correlation_scale)
    y = np.arange(-correlation_scale, correlation_scale)
    X, Y = np.meshgrid(x, y)
    dist = np.sqrt(X*X + Y*Y)
    filter_kernel = np.exp(-dist**2/(2*correlation_scale))
    noise = np.random.normal(0.0, 1.0, size=(ny,nx)) #.reshape(ny,nx)
    noise = scipy.signal.fftconvolve(noise, filter_kernel, mode='same')
    noise = np.real(noise)*(1.0/(np.std(np.real(noise))+1e-20))
    # noise = np.fft.fft2(noise)
    # filter_kernel = np.fft.fft2(filter_kernel)
    # noise = np.fft.ifft2(noise*filter_kernel)
    # noise = np.real(noise)*(1.0/(np.std(np.real(noise))+1e-20))
    return noise
######################################
def spatially_correlated_random_init(num,nx,ny,tau=50.0):
    """
    Spatially correlated random variable 
    sptial demensions are [nx,ny]
    tau_s = decorrelation length
    """
    np.random.seed(num)
    correlation_scale = tau
    x = np.arange(-correlation_scale, correlation_scale)
    y = np.arange(-correlation_scale, correlation_scale)
    X, Y = np.meshgrid(x, y)
    dist = np.sqrt(X*X + Y*Y)
    filter_kernel = np.exp(-dist**2/(2*correlation_scale))
    noise = np.random.normal(0.0, 1.0, size=(ny,nx)) #.reshape(ny,nx)
    noise = scipy.signal.fftconvolve(noise, filter_kernel, mode='same')
    noise = np.real(noise)*(1.0/(np.std(np.real(noise))+1e-20))
    # noise = np.fft.fft2(noise)
    # filter_kernel = np.fft.fft2(filter_kernel)
    # noise = np.fft.ifft2(noise*filter_kernel)
    # noise = np.real(noise)*(1/(np.std(np.real(noise))+1e-20))
    return noise
######################################
def make_runoff(inputlist):
    """
    Create new runoff file 
    """
    # print (inputlist)
    yyyy     =inputlist[0]
    mm       =inputlist[1]
    dd       =inputlist[2]
    ens_char =inputlist[3]
    coeff    =float(inputlist[4])
    std      =float(inputlist[5])
    ifile    =inputlist[6]
    roff_mean=inputlist[7]
    nx       =int(inputlist[8])
    ny       =int(inputlist[9])
    ofile    =inputlist[10]
    ifile=np.fromfile(ifile,np.float32).reshape(ny,nx)
    roff_mean=np.fromfile(roff_mean,np.float32).reshape(ny,nx)
    # roff_total=np.fromfile(outdir+"/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
    tmp=ifile*coeff + roff_mean*std
    tmp.astype(np.float32)
    tmp.tofile(ofile)
    return 0
######################################
def make_runoff_temporal(inputlist):
    """
    make spatially and temporally correlated runoff: lognormal distribution
    Refferce:
    Nijssen, B., & Lettenmaier, D. P. (2004). Effect of precipitation 
    sampling error on simulated hydrological fluxes and states: Anticipating 
    the Global Precipitation Measurement satellites. Journal of Geophysical 
    Research: Atmospheres, 109(2), 1–15. https://doi.org/10.1029/2003jd003497
    """
    syear=int(inputlist[0])
    eyear=int(inputlist[1])
    ens_char=inputlist[2]
    runname=inputlist[3]
    alpha=float(inputlist[4])
    beta=float(inputlist[5])
    E=float(inputlist[6])
    tau_s=float(inputlist[7])
    rundir=inputlist[8]
    outdir=inputlist[9]
    # get runoff metadata
    nx,ny,prefix,suffix=get_runoff_metadata(runname)
    # create necessary varibales
    # get spatially correlated normal distributiion
    num=int(ens_char)
    s=spatially_correlated_random_init(num,nx,ny,tau_s)#*rand
    # s=s*(1/(np.std(s)+1e-20))
    # date
    start_dt=datetime.date(syear,1,1)
    end_dt=datetime.date(eyear,12,31)
    start=0
    last=(end_dt-start_dt).days + 1
    for day in np.arange(start,last):
        target_dt=start_dt+datetime.timedelta(days=day)
        yyyy='%04d' % (target_dt.year)
        mm='%02d' % (target_dt.month)
        dd='%02d' % (target_dt.day)
        # get spatially correlated normal distributiion
        # w=spatially_correlated_random(nx,ny)
        # stochastic term with mean 0 and variance 1 following a Gaussian distribution
        w=np.random.normal(0,1)
        # print ("w ", np.max(w),np.min(w))
        s=alpha*s+np.sqrt(1-alpha**2)*w
        #===
        Rcoff = ((1.0 + beta) / np.sqrt(E**2 + 1.0)) * np.exp(np.sqrt(np.log(E**2 + 1)) * s)
        Rcoff = np.float32(Rcoff)
        # Rcoff.astype(np.float32, casting='unsafe', copy=True)
        # Rcoff=np.zeros([ny,nx],np.float32)
        # for ix in np.arange(nx):
        #     for iy in np.arange(ny):
        #         Rcoff[iy,ix]=((1.0+beta)/math.sqrt(E**2+1.0))*math.exp(math.sqrt(E**2+1)*s[iy,ix])
        # ifile=rundir+"/"+prefix+yyyy+mm+dd+suffix
        ifile=np.fromfile(rundir+"/"+prefix+yyyy+mm+dd+suffix,np.float32).reshape(ny,nx)
        ifile=np.nan_to_num(ifile)
        ifile=np.float32(ifile)
        # ifile.astype(np.float32, casting='unsafe', copy=True)
        # Rc=np.zeros([ny,nx],np.float32)
        Rc=Rcoff*ifile
        # Rc=Rc*(Rc>1e-2)*1.0
        # Rc.astype(np.float32, casting='unsafe', copy=True)
        Rc=np.abs(Rc) #(Rc>0.0)*1.0 #remove negative values
        Rc=np.float32(Rc)
        # print (yyyy,mm,dd,ens_char,np.max(Rc),np.min(Rc)) #np.max(Rcoff),np.min(Rcoff),np.max(ifile),np.min(ifile),
        ofile=outdir+"/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
        Rc.tofile(ofile)
    return 0
######################################
def make_runoff_normal(inputlist):
    """
    make spatially and temporally correlated runoff: normal distribution
    """
    syear=int(inputlist[0])
    eyear=int(inputlist[1])
    ens_char=inputlist[2]
    runname=inputlist[3]
    alpha=float(inputlist[4])
    beta=float(inputlist[5])
    E=float(inputlist[6])
    rand=float(inputlist[7])
    rundir=inputlist[8]
    outdir=inputlist[9]
    # get runoff metadata
    nx,ny,prefix,suffix=get_runoff_metadata(runname)
    # create necessary varibales
    num=int(ens_char)
    s=spatially_correlated_random_init(num,nx,ny)#*rand
    # s=s*(1/(np.std(s)+1e-20))
    # date 
    start_dt=datetime.date(syear,1,1)
    end_dt=datetime.date(eyear,12,31)
    start=0
    last=(end_dt-start_dt).days + 1
    for day in np.arange(start,last):
        target_dt=start_dt+datetime.timedelta(days=day)
        yyyy='%04d' % (target_dt.year)
        mm='%02d' % (target_dt.month)
        dd='%02d' % (target_dt.day)
        # get spatially correlated normal distributiion
        w=spatially_correlated_random(nx,ny)
        # stochastic term with mean 0 and variance 1 following a Gaussian distribution
        # w=np.random.normal(0,1)
        # print ("w ", np.max(w),np.min(w))
        s=alpha*s+np.sqrt(1-alpha**2)*w
        #===
        Rcoff = (1.0 + s) #((1.0 + beta) / np.sqrt(E**2 + 1.0)) * np.exp(np.sqrt(np.log(E**2 + 1)) * s)
        Rcoff = np.float32(Rcoff)
        # Rcoff.astype(np.float32, casting='unsafe', copy=True)
        # Rcoff=np.zeros([ny,nx],np.float32)
        # for ix in np.arange(nx):
        #     for iy in np.arange(ny):
        #         Rcoff[iy,ix]=((1.0+beta)/math.sqrt(E**2+1.0))*math.exp(math.sqrt(E**2+1)*s[iy,ix])
        # ifile=rundir+"/"+prefix+yyyy+mm+dd+suffix
        ifile=np.fromfile(rundir+"/"+prefix+yyyy+mm+dd+suffix,np.float32).reshape(ny,nx)
        ifile=np.nan_to_num(ifile)
        ifile=np.float32(ifile)
        # ifile.astype(np.float32, casting='unsafe', copy=True)
        # Rc=np.zeros([ny,nx],np.float32)
        Rc=Rcoff*ifile
        # Rc=Rc*(Rc>1e-2)*1.0
        # Rc.astype(np.float32, casting='unsafe', copy=True)
        Rc=np.abs(Rc) #(Rc>0.0)*1.0 #remove negative values
        Rc=np.float32(Rc)
        # print (yyyy,mm,dd,ens_char,np.max(Rc),np.min(Rc)) #np.max(Rcoff),np.min(Rcoff),np.max(ifile),np.min(ifile),
        ofile=outdir+"/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
        Rc.tofile(ofile)
    return 0
######################################
def copy_runoff_subbasin(inputlist): #do it parallel
    iname=inputlist[0]
    oname=inputlist[1]
    ens  =int(inputlist[2]) 
    nx   =int(inputlist[3])
    ny   =int(inputlist[4])
    # distopen=float(inputlist[2])
    runoff=np.fromfile(iname,np.float32).reshape(ny,nx)
    runoff=(np.ma.masked_equal(runoff,-9999.0)*rand).filled(-9999.0)
    runoff.tofile(oname)
    return 0
######################################
def corr_norm(sigma,n_samples,dcl=500):
    # make spatially correlated normal random fields
    # find covariances
    # l=40
    # decorrelation length
    dcl=500 #km
    sigma=0.01
    # CaMa-Flood lat lon
    lonlat=pm.CaMa_dir()+"/map/"+pm.mapname()+"/lonlat.bin"
    lonlat=np.fromfile(lonlat,np.float32).reshape(2,250,350)
    fname="outsub_amz_06min.txt"
    # fname="outsub.txt"
    ids=[]
    lxx=[]
    lyy=[]
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line = filter(None,re.split(" ", line))
        ids.append(float(line[0]))
        lxx.append(int(line[1]))
        lyy.append(int(line[2]))
    #---------------------------
    l=len(ids)
    # print "/src/make_covariance "+str(l)+" "+str(1)+" "+pm.CaMa_dir()+"/ "+fname
    # os.system("./src/make_covariance "+str(l)+" "+str(1)+" "+pm.CaMa_dir()+"/ "+fname)
    # calculate covariance
    print ("calculate covariance", l, "x", l)
    mean=np.ones([l],np.float32)*sigma #*pm.corruptman_base()
    covariance=np.zeros([l,l],np.float32)
    # covariance=np.fromfile("cov.bin",np.float32).reshape(l,l)
    for i in np.arange(l):
        for j in np.arange(l):
            ix1 = lxx[i] - 1
            iy1 = lyy[i] - 1
            ix2 = lxx[j] - 1
            iy2 = lyy[j] - 1
            lon1=lonlat[0,iy1,ix1]
            lat1=lonlat[1,iy1,ix1]
            lon2=lonlat[0,iy2,ix2]
            lat2=lonlat[1,iy2,ix2]
            lag=hubeny(lat1,lon1,lat2,lon2)*1e-3
            covariance[i,j]=covariance_lag(lag,sigma,dcl)
            # print (lag, covariance[i,j])
    #--
    # # Compute the Cholesky decomposition
    # c=spla.cholesky(covariance, lower=True)
    # n_samples=l
    # create multivarite samples 
    L=spla.cholesky(covariance)
    Z=np.random.normal(size=(n_samples,covariance.shape[0]))
    # print Z.dot(L) + mean
    # print np.random.multivariate_normal(mean,covariance,size=1) + mean
    return np.abs(Z.dot(L) + mean)
    # return np.random.multivariate_normal(mean,covariance,size=1) + mean
    # return np.random.default_rng().multivariate_normal(mean,covariance,size=1,method='eigh') + mean
    # #--
    # return np.dot(c,mean)
######################################
def hubeny(lat1,lon1,lat2,lon2):
    pi = math.pi #math.atan(1.0)*4.0
    a  = 6378137
    b  = 6356752.314140
    e2 = 0.00669438002301188
    a_1_e2 = 6335439.32708317
    latrad1   = lat1 * pi / 180.0
    latrad2   = lat2 * pi / 180.0
    lonrad1   = lon1 * pi / 180.0
    lonrad2   = lon2 * pi / 180.0
    #-----------------------------------
    latave    = (latrad1 + latrad2)/2.0
    dlat      = latrad2 - latrad1
    dlon      = lonrad2 - lonrad1
    #-----------------------------------
    dlondeg   = lon2 - lon1
    if abs(dlondeg) > 180.0:
        dlondeg = 180.0 - abs(dlondeg)%180.0
        dlon    = dlondeg * pi / 180.0
    #-----------------------------------
    W  = math.sqrt(1.0 - e2 * math.sin(latave)**2.0 )
    M  =  a_1_e2 / (W**3.0)
    N  =  a / W
    hubeny_real  = math.sqrt( (dlat * M)**2.0 + (dlon * N * math.cos(latave))**2.0 )
    return hubeny_real
######################################
def covariance_lag(lag,sigma,L):
    return (sigma**2)*(math.exp(-1.0*lag/L))
######################################
# def spatiotemp_random():
#     #--
#     sigma=1.0
#     n_samples=20
#     dcl=500 #km
#     dt=1.0
#     tau=10.0
#     rho=temp_para(dt,tau)
#     w=corr_norm(sigma,n_samples,dcl)
#     s=np.random.normal(0.0,1.0,size=(n_samples,np.shape(w)[1]))
#     st=rho*s*+math.sqrt(1-rho)*w
#     u=0.5*math.erfc(st/math.sqrt(2))
######################################
def temp_para(dt,tau):
    return 1-(dt/(tau+1.0e-20))
######################################
def get_runoff_metadata(runname):
    if runname=="ELSE_KIM2009":
        nx,ny=360,180
        prefix="Roff____"
        suffix=".one"
    if runname=="E2O":
        nx,ny=1440,720
        prefix="Roff__"
        suffix=".one"
    if runname=="ERA20CM":
        nx,ny=360,180
        prefix="Roff__"
        suffix=".one"
    if runname=="isimip3a":
        nx,ny=720,360
        prefix="G5oFNRCD"
        suffix=".hlf"
    if runname=="ERA5":
        nx,ny=3600,1800
        prefix="Roff____"
        suffix=".sixmin"
    return nx,ny,prefix,suffix
######################################
def val_distopen():
    return 1.0
######################################
def val_diststd():
    return 0.25
######################################
def val_tau_t():
    return 150 # decorrelation time
######################################
def val_tau_s():
    return 50 # spatial decorrelation length
######################################
def val_alpha():
    return 1.0 - (1.0/(val_tau_t()+1e-20))
######################################
def val_beta():
    return 0.0
######################################
def val_E():
    return 0.30
######################################
def runoff_name(): 
    # return "isimip3a"
    return "ERA5"
######################################
def runoff_dir():
    # return "/work/a04/julien/CaMa-Flood_v4/inp/isimip3a/runoff"
    return "/work/a02/menaka/ERA5/bin"
######################################
def out_dir():
    return "/work/a06/menaka/ensemble_simulations/CaMa_in"
######################################
if __name__ == "__main__":
    # prepare_input()
    # print (corr_norm(1.00,20))
    sys.path.append("../gosh/")
    # import params_virt as pm
    import params_real as pm
    # simple
    prep_runoff_ensemble(val_distopen(),val_diststd(),20,runoff_name(),runoff_dir(),out_dir(),2000,2010,"simple",beta=val_beta(),E=val_E(),aplha=val_alpha())
    # lognormal
    prep_runoff_ensemble(val_distopen(),val_diststd(),20,runoff_name(),runoff_dir(),out_dir(),2000,2010,"lognormal",beta=val_beta(),E=val_E(),aplha=val_alpha())
    # normal
    prep_runoff_ensemble(val_distopen(),val_diststd(),20,runoff_name(),runoff_dir(),out_dir(),2000,2010,"normal",beta=val_beta(),E=val_E(),aplha=val_alpha())
    