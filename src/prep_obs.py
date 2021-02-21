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
def mk_dir(sdir):
    try:
        os.makedirs(sdir)
    except:
        pass
#########################
def slink(src,dst):
  try:
    os.symlink(src,dst)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST:
      os.remove(dst)
      os.symlink(src,dst)
    else:
      raise
#########################
def get_HydroWeb():
	#---
	lname=[]
	xlist=[]
	ylist=[]
	lEGM08=[]
	lEGM96=[]
	satellite=[]
	fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+pm.mapname()+".txt"
	with open(fname,"r") as f:
		lines=f.readlines()
	for line in lines[1::]:
		line    = re.split(" ",line)
		line    = list(filter(None, line))
		#print line
		num     = line[0]
		station = line[1]
		riv     = re.split("_",station)[1]
		lon     = float(line[2])
		lat     = float(line[3])
		ix      = int(line[4])
		iy      = int(line[5])
		EGM08   = float(line[8])
		EGM96   = float(line[9])
		sat     = line[10].split()[0]
		#------
		lname.append(station)
		xlist.append(ix)
		ylist.append(iy)
		lEGM08.append(EGM08)
		lEGM96.append(EGM96)
		satellite.append(sat)
	return lname, xlist, ylist, lEGM08, lEGM96, satellite
#########################
def write_txt(inputlist):
	yyyy=inputlist[0]
	mm=inputlist[1]
	dd=inputlist[2]
	HydroWeb_dir=inputlist[3]
	target_dt=datetime.date(int(yyyy),int(mm),int(dd))
	txtfile=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/"+yyyy+mm+dd+".txt"
	# print txtfile
	pnum=len(lname)
	# print pnum
	with open(txtfile,"w") as txtf:
		for point in np.arange(pnum):
			#--
			#print name
			lwse=[]
			wseo=-9999.0
			#--read HydroWeb data
			iname=HydroWeb_dir+"/data/hydroprd_"+lname[point]+".txt"
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
				wse  = float(ll_hyd[2]) + lEGM08[point] - lEGM96[point]
				lwse.append(wse)
				now  = datetime.date(year,mon,da)
				if now == target_dt:
					wseo=wse
			if wseo == -9999.0:
				continue
			# write txt file
			iix=xlist[point]
			iiy=ylist[point]
			mean_wse=np.mean(np.array(lwse))
			std_wse=np.std(np.array(lwse))
			sat=satellite[point]
			line="%04d	%04d	%10.4f	%10.4f	%10.4f	%s\n"%(iix,iiy,wseo,mean_wse,std_wse,sat)
			txtf.write(line)
			# print (line)
			return 0
#########################
def prepare_obs():
	#
	global lname, xlist, ylist, lEGM08, lEGM96, satellite
	lname, xlist, ylist, lEGM08, lEGM96, satellite = get_HydroWeb()
	syear,smon,sday=pm.starttime()
	eyear,emon,eday=pm.endtime()
	start_dt=datetime.date(syear,smon,eday)
	end_dt=datetime.date(eyear,emon,eday)
	obs_dir="/cluster/data6/menaka/HydroWeb"
	if pm.obs_name() == "HydroWeb":
		obs_dir="/cluster/data6/menaka/HydroWeb"
	start=0
	last=(end_dt-start_dt).days
	#-------
	inputlist=[]
	for day in np.arange(start,last):
		target_dt=start_dt+datetime.timedelta(days=day)
		yyyy='%04d' % (target_dt.year)
		mm='%02d' % (target_dt.month)
		dd='%02d' % (target_dt.day)
		print yyyy,mm,dd,obs_dir
		inputlist.append([yyyy,mm,dd,obs_dir])
	# write text files parallel
	p=Pool(20)
	p.map(write_txt,inputlist)
	p.terminate()
	return 0
####################################
# if __name__ == "__main__":
#     prepare_obs()