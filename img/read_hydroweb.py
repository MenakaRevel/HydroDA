#!/usr/bin/python
import os
import re
import datetime
import numpy as np
################################
# Fuctions to read HydroWeb data
################################
def hydroweb_river_name(mapname="glb_15min"):
    # directory
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    river=[]
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        print line
        station = line[1]
        riv     = re.split("_",station)[1]
        river.append(riv)
    return river
################################
def get_hydroweb(mapname="glb_15min"):
    # directory
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    river=[]
    pname=[]
    xlist=[]
    ylist=[]
    egm_d=[]
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        print line
        station = line[1]
        riv     = re.split("_",station)[1]
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        river.append(riv)
        pname.append(station)
        xlist.append(ix)
        ylist.append(iy)
        egm_d.append(EGM96-EGM08)
    return river,pname,xlist,ylist,egm_d
##################################
def get_hydroweb_loc(rivername,mapname="glb_15min"):
    # directory
    fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    #river=[]
    pname=[]
    xlist=[]
    ylist=[]
    egm08=[]
    egm96=[]
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        #print line
        station = line[1]
        riv     = re.split("_",station)[1]
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        if rivername==riv:
            #river.append(riv)
            pname.append(station)
            xlist.append(ix)
            ylist.append(iy)
            egm08.append(EGM08)
            egm96.append(EGM96)
    return pname,xlist,ylist,egm08,egm96
##################################
def HydroWeb_WSE(station,syear,eyear,smon=1,emon=12,sday=1,eday=31):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read hydroweb
    #station="R_con_con_env_0429_01"
    #satellite=station.split("_")[2]
    #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
    fname="/cluster/data6/menaka/HydroWeb/data/hydroprd_"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=33
    #--
    time=[] # time in days
    data=[] # WSE in [m]
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = re.split(" ",line)
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        wse  = float(line[2])
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        if now < start and now > end:
            continue
        data.append(wse)
        lag  = int((now-start).days)
        time.append(lag)
    return time, data
#####################################
def HydroWeb_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read hydroweb
    #station="R_con_con_env_0429_01"
    #satellite=station.split("_")[2]
    #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
    fname="/cluster/data6/menaka/HydroWeb/data/hydroprd_"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=33
    #--
    time=int((end-start).days + 1) # time in days
    data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = re.split(" ",line)
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        wse  = float(line[2])
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        #print yyyy, mm, dd, int((now-start).days), wse
        if now > start and now < end:
            nowtime=int((now-start).days)
            data[nowtime]=wse+egm08-egm96
    return data
#####################################
def altimetry(name,mapname="glb_15min"):
    # directory
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        station = line[1]
        if station==name:
            alti= float(line[6])
    return alti
################################
