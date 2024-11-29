#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
import json
import pandas as pd
################################
# Functions to read CGLS data
################################
def get_cgls_locs(mapname="glb_15min"):
    # directory
    fname="/cluster/data6/menaka/CGLS/data/river/CGLS_alloc_"+mapname+".txt"
    #--
    nums=[]
    river=[]
    pname=[]
    lons =[]
    lats =[]
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
        num     = line[0]
        station = line[1]
        riv     = re.split("_",station)[1]
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        #-----------------------
        nums.append(num)
        river.append(riv)
        pname.append(station)
        lons.append(lon)
        lats.append(lat)
        xlist.append(ix)
        ylist.append(iy)
        egm08.append(EGM08)
        egm96.append(EGM96)
    return nums,river,pname,lons,lats,xlist,ylist,egm08,egm96
################################
def get_cgls_loc(rivername,mapname="conus_06min",fname="/cluster/data6/menaka/HydroDA/dat/CGLS_alloc_conus_06min_org.txt"):
    # read metadata
    meta="/work/a06/menaka/CGLS/data/river/CGLS_ALTI_V2.1.0.txt"
    dfmeta=pd.read_csv(meta,sep=";")
    names=dfmeta[dfmeta["Basin"]==rivername]["Station"].values
    #read fname
    dfcgls=pd.read_csv(fname,delim_whitespace=True)
    pname=dfcgls[dfcgls["station"].isin(names)]["station"].values
    xlist=dfcgls[dfcgls["station"].isin(names)]["ix"].values - 1
    ylist=dfcgls[dfcgls["station"].isin(names)]["iy"].values - 1
    egm08=dfcgls[dfcgls["station"].isin(names)]["EGM08"].values
    egm96=dfcgls[dfcgls["station"].isin(names)]["EGM96"].values
    return pname,xlist,ylist,egm08,egm96
################################
def cgls_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    #---
    time=int((end-start).days + 1) # time in days
    data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
    # read CGLS
    fname="/work/a06/menaka/CGLS/data/river/"+station+".json"
    with open(fname) as f:
        alldata    = json.load(f)
        cgls_data  = alldata["data"]
    #----------------------------
    for line in range(len(cgls_data)):
        date    = cgls_data[line]["datetime"]
        date    = re.split(" ",date)[0]
        date    = re.split("/",date)
        yyyy    = int(date[0])
        mm      = int(date[1])
        dd      = int(date[2])
        now     = datetime.date(yyyy,mm,dd)
        wse     = cgls_data[line]["water_surface_height_above_reference_datum"]
        if now > start and now < end:
            nowtime=int((now-start).days)
            data[nowtime]=wse+egm08-egm96
    return data
#####################################
def cgls_WSE(station,syear,eyear,smon=1,emon=12,sday=1,eday=31):
    start = datetime.date(syear, smon, sday)
    end = datetime.date(eyear, emon, eday)
    # time = (end - start).days + 1  # time in days
    # data = np.full(time, -9999.0, dtype=np.float32)  # WSE in [m], -9999.0 for no observations

    fname = "/work/a06/menaka/CGLS/data/river/" + station + ".json"
    with open(fname) as f:
        alldata = json.load(f)
        cgls_data = alldata["data"]

    date_format = "%Y/%m/%d"
    # start_time = (start - datetime.date(2000, 1, 1)).days  # Start time relative to a reference date

    # Create a DataFrame from cgls_data
    df = pd.DataFrame(cgls_data)
    df["datetime"] = pd.to_datetime(df["datetime"].str.split(" ").str[0], format=date_format)

    # Filter and update data efficiently using DataFrame operations
    mask = (df["datetime"] >= start) & (df["datetime"] <= end)
    df_filtered = df[mask].copy()
    df_filtered["days_from_start"] = (df_filtered["datetime"] - start).dt.days

    # Separate data and time into arrays
    time_array = df_filtered["days_from_start"].values
    data_array = df_filtered["water_surface_height_above_reference_datum"].values #+ egm08 - egm96

    return time_array, data_array
#####################################
def metadata():
    # directory 
    hydroweb="/cluster/data6/menaka/CGLS/data/river/CGLS_Station_list.txt"
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    sta={}
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        id      = int(line[0])
        station = line[1]
        river   = line[2]
        basin   = line[3]
        country = line[4]
        sat     = line[10]
        startdt = line[11]
        enddt   = line[13]
        status  = line[15].split("\n")[0]
        #----
        sta[station]=[river,basin,country,sat,startdt,enddt,status]
    #----
    return sta
################################