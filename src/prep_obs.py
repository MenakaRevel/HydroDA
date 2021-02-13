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
def get_GRDC():
    #---
    lname=[]
    rivnm=[]
    xlist=[]
    ylist=[]
    fname=pm.CaMa_dir()+"/map/"+pm.mapname()+"/grdc_loc.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print line
        num     = line[0].strip()
        river   = line[1].strip()
        station = line[2].strip()
        ix      = int(line[3])
        iy      = int(line[4])
        if int(line[5])!=-9999:
            continue
        #------
        print (river, station)
        lname.append(station)
        rivnm.append(river)
        xlist.append(ix)
        ylist.append(iy)
    return lname, rivnm, xlist, ylist
#########################
def prepare_obs():
    #### define globally
    global sim, shared_array_sim
    global start_dt,end_dt
    global lname, rivnm, xlist, ylist
    global start, last
    global means, stds
    global pnum
    global nx, ny, gsize
    nx, ny, gsize = pm.map_dimension()
    # print (nx,ny,gsize)
    syear,smon,sday=pm.starttime()
    eyear,emon,eday=pm.endtime()
    start_dt=datetime.date(syear,smon,eday)
    end_dt=datetime.date(eyear,emon,eday)
    start=0
    last=(end_dt-start_dt).days
    #-------
    lname, rivnm, xlist, ylist=get_GRDC()
    pnum=len(lname)
    print (pnum)
    # multiprocessing array
    sim=np.ctypeslib.as_ctypes(np.zeros([last,pnum],np.float32))
    shared_array_sim  = sharedctypes.RawArray(sim._type_, sim)
    
    inputlist=[]
    chpnum='%d'%(pnum)
    for year in np.arange(syear,eyear):
        yyyy='%04d' % (year)
        inputlist.append([yyyy]) #,chpnum
    #---------------------
    # map(read_data, inputlist)
    p   = Pool(20)
    res = p.map(read_data, inputlist)
    sim = np.ctypeslib.as_array(shared_array_sim)
    p.terminate()
    #---------------------
    means=np.mean(sim,axis=0)
    stds=np.std(sim,axis=0)
    #-------
    inputlist=[]
    for day in np.arange(start,last+1):
        target_dt=start_dt+datetime.timedelta(days=day)
        yyyy='%04d' % (target_dt.year)
        mm='%02d' % (target_dt.month)
        dd='%02d' % (target_dt.day)
        chday= '%02d' %(day)
        chpnum='%d'  % (pnum)
        inputlist.append([yyyy,mm,dd,chday]) #,chpnum
    #---------------------
    # map(write_text, inputlist)
    p = Pool(20)
    p.map(write_text, inputlist)
    p.terminate()
###########################
def write_text(inputlist):
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    day=int(inputlist[3])
    # pnum=int(inputlist[4])
    txtfile="./assim_out/obs/"+yyyy+mm+dd+".txt"
    print (txtfile)
    with open(txtfile,"w") as txtf:
        for point in np.arange(pnum):
            #### ix    iy    value    mean    std
            line="%04d    %04d    %8.4f    %8.4f    %8.4f\n"%(xlist[point],ylist[point],sim[day,point],means[point],stds[point])
            txtf.write(line)
            print (line)
#########################
def get_HydroWeb(inputlist):
	yyyy=inputlist[0]
	mm=inputlist[1]
	dd=inputlist[2]
	obs_dir="/cluster/data6/menaka/ensemble_org/CaMa_out/AMZE2O003" # folder where Simulated discharge
	target_dt=datetime.date(int(yyyy),int(mm),int(dd))
	txtfile="./assim_out/obs/"+yyyy+mm+dd+".txt"
	print txtfile
	with open(txtfile,"w") as txtf:
		for name in lname:
			#--
			#print name
			lwse=[]
			wseo=-9999.0
			#--read HydroWeb data
			iname=HydroWeb_dir+"/data/hydroprd_"+name+".txt"
			with open(iname,"r") as f_hyd:
				l_hyd=f_hyd.readlines()
			for ll_hyd in l_hyd[33::]:
				ll_hyd = re.split(" ",ll_hyd)
				ll_hyd = list(filter(None, ll_hyd))
				date = ll_hyd[0]
				date = re.split("-",date)
				year = int(date[0])
				mon  = int(date[1])
				da   = int(date[2])
				wse  = float(ll_hyd[2]) + lEGM08[name] - lEGM96[name]
				lwse.append(wse)
				now  = datetime.date(year,mon,da)
				if now == target_dt:
					wseo=wse
			if wseo == -9999.0:
				continue
			# write txt file
			iix=xlist[name]
			iiy=ylist[name]
			mean_wse=np.mean(np.array(lwse))
			std_wse=np.std(np.array(lwse))
			sat=satellite[name]
			line="%04d	%04d	%8.4f	%8.4f	%8.4f	%s\n"%(iix,iiy,wseo,mean_wse,std_wse,sat)
			txtf.write(line)
			print name, line
####################################
def read_data(inputlist):
    yyyy = inputlist[0]
    #pnum = int(inputlist[1])
    #odir = inputlist[1]
    print (yyyy)
    #--
    tmp_sim  = np.ctypeslib.as_array(shared_array_sim)

    # year, mon, day
    year=int(yyyy)
    
    if calendar.isleap(year):
        dt=366
    else:
        dt=365

    # timimgs
    target_dt=datetime.date(year,1,1)
    st=(target_dt-start_dt).days
    et=st+dt
    if et >= last:
        et=None

    # simulated discharge
    # nx, ny = 1440,720
    print (nx,ny)
    fname=pm.obs_dir()+"/outflw"+yyyy+".bin"
    # fname="/cluster/data6/menaka/ensemble_org/CaMa_out/GLBE2O003/outflw"+yyyy+".bin"
    simfile=np.fromfile(fname,np.float32).reshape([dt,ny,nx])
    #-------------
    for point in np.arange(pnum):
        ix,iy=xlist[point],ylist[point]
        tmp_sim[st:et,point]=simfile[:,iy-1,ix-1]
        # ix1,iy1,ix2,iy2=x1list[point],y1list[point],x2list[point],y2list[point]
        # if ix2 == -9999 or iy2 == -9999:
        #     tmp_sim[st:et,point]=simfile[:,iy1-1,ix1-1]
        # else:
        #     tmp_sim[st:et,point]=simfile[:,iy1-1,ix1-1]+simfile[:,iy2-1,ix2-1]
####################################
if __name__ == "__main__":
    prepare_obs()