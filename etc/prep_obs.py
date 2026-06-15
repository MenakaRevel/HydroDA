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
import sys

# #external python codes
# dir_param="../gosh"
# sys.path.append(dir_param)
# import params as pm
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
	#=============
	lname=[]
	xlist=[]
	ylist=[]
	leledf=[]
	lEGM08=[]
	lEGM96=[]
	satellite=[]
	# fname=pm.DA_dir()+"/dat/HydroWeb_alloc_"+pm.mapname()+".txt"
	# fname=pm.DA_dir()+"/dat/HydroWeb_alloc_"+pm.mapname()+"_new.txt"
	# fname=pm.DA_dir()+"/dat/HydroWeb_alloc_"+pm.mapname()+"_amz.txt"
	# fname=obs_list()
	fname=HydroWeb_list()
	#=========================
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
		eledif  = float(line[7])
		EGM08   = float(line[8])
		EGM96   = float(line[9])
		sat     = line[10].split()[0]
		# print (riv,station,sat)
		#------
		lname.append(station)
		xlist.append(ix)
		ylist.append(iy)
		leledf.append(eledif)
		lEGM08.append(EGM08)
		lEGM96.append(EGM96)
		satellite.append(sat)
	return lname, xlist, ylist, leledf, lEGM08, lEGM96, satellite
#########################
def read_HydroWeb(yyyy,mm,dd,name,EGM08,EGM96,eledf=0.0):
	target_dt=datetime.date(int(yyyy),int(mm),int(dd))
	HydroWeb_dir="/cluster/data6/menaka/HydroWeb"
	#print name
	lwse=[]
	wseo=-9999.0
	#--read HydroWeb data
	iname=HydroWeb_dir+"/data/hydroprd_"+name+".txt"
	# print (iname)
	with open(iname,"r") as f_hyd:
		l_hyd=f_hyd.readlines()
	for ll_hyd in l_hyd[33::]:
		ll_hyd = re.split(" ",ll_hyd)
		ll_hyd = list(filter(None, ll_hyd))
		date = ll_hyd[0]
		date = re.split("-",date)
		year = int(date[0])
		mon  = int(date[1])
		day  = int(date[2])
		wse  = float(ll_hyd[2]) + EGM08 - EGM96 + eledf
		lwse.append(wse)
		now  = datetime.date(year,mon,day)
		# print (yyyy, mm, dd)
		# print ("data:", year,mon,day,wse)
		if now == target_dt:
			wseo=wse
	if wseo==-9999.0:
		mean_wse=-9999.0
		std_wse=-9999.0
	else:
		mean_wse=np.mean(np.array(lwse))
		std_wse=np.std(np.array(lwse))
	# print (yyyy, mm, dd, wseo, mean_wse, std_wse)
	return wseo, mean_wse, std_wse
#########################
def HydroWeb_data(yyyy,mm,dd):
	lname =[]
	xlist =[]
	ylist =[]
	l_wse =[]
	m_wse =[]
	s_wse =[]
	l_sat =[]
	lstan, xcods, ycods, leledif, lEGM08, lEGM96, satellite = get_HydroWeb()
	pnum=len(lstan)
	# print (pnum)
	for point in np.arange(pnum):
		# # == read relevant observation data ==
		# if pm.obs_name() == "HydroWeb":
		# 	# == for HydroWeb data ==
		wseo, mean_wse, std_wse = read_HydroWeb(yyyy,mm,dd,lstan[point],lEGM08[point],lEGM96[point])#,leledif[point])
		# else:
		# 	wseo, mean_wse, std_wse = read_HydroWeb(yyyy,mm,dd,lname[point],lEGM08[point],lEGM96[point],leledif[point])
		# print (point, wseo)
		if wseo == -9999.0:
			continue
		# print (point, wseo)
		iix=xcods[point]
		iiy=ycods[point]
		sat=satellite[point]
		#====================
		xlist.append(iix)
		ylist.append(iiy)
		l_wse.append(wseo)
		m_wse.append(mean_wse)
		s_wse.append(std_wse)
		l_sat.append(sat)
	return xlist, ylist, l_wse, m_wse, s_wse, l_sat
####################################
def swot_data(yyyy,mm,dd):
	# prepare sythetic observations using
	# pre-simulated data
	# river width thershold
	rivwdth_thr=50.0 #m
	nx,ny,gsize = map_dimension()
	ny_swot = min(ny,640)
	day=SWOT_day(yyyy,mm,dd)
	SWOTDD="%02d"%(day)
	fname="../sat/mesh_day"+SWOTDD+".bin" # for glb_15min
	mesh_in=np.fromfile(fname,np.float32).reshape([ny_swot,nx])
	mesh=(mesh_in>=10)*(mesh_in<=60)
	meshP=mesh-1000*(mesh<0.1)
	SWOTmesh=np.zeros([ny,nx],np.float32)
	SWOTmesh[40:680,:]=meshP
	fname=CaMa_dir()+"/map/"+mapname()+"/rivwth_gwdlr.bin"
	rivwth=np.fromfile(fname,np.float32).reshape(ny,nx)
	damloc=dam_loc()
	obs=(SWOTmesh>=1.0)*(rivwth>=rivwdth_thr)*(damloc>0.0)*1.0
	lname =[]
	xlist =[]
	ylist =[]
	l_wse =[]
	m_wse =[]
	s_wse =[]
	l_sat =[]
	# leledif, lEGM08, lEGM96, satellite
	#===================================
	# odir="/cluster/data6/menaka/ensemble_org/CaMa_out/E2O003"
	odir=obs_dir()
	fname=odir+"/sfcelv"+yyyy+".bin"
	year=int(yyyy)
	mon=int(mm)
	day=int(dd)
	if calendar.isleap(year):
		nt=366
	else:
		nt=365
	orgfile=np.fromfile(fname,np.float32).reshape([nt,ny,nx])
	obs_err=SWOT_observation_error()
	#-----------------------
	start_dt=datetime.date(year,1,1)
	target_dt=datetime.date(year,mon,day)
	it=(target_dt-start_dt).days
	for ix in np.arange(nx):
		for iy in np.arange(ny):
			if obs[iy,ix] == 1.0:
				wse=orgfile[it,iy,ix] + err_rand(obs_err,ix,iy)
				# print (ix,iy,wse[0])
				l_wse.append(wse[0])
				xlist.append(ix+1)
				ylist.append(iy+1)
				m_wse.append(-9999.0)
				s_wse.append(-9999.0)
				l_sat.append("SWOT")
	return xlist, ylist, l_wse, m_wse, s_wse, l_sat
####################################
def SWOT_day(yyyy,mm,dd):
	st_year,st_month,st_date=starttime()
	start_time=datetime.date(st_year,st_month,st_date)
	this_time=datetime.date(int(yyyy),int(mm),int(dd))
	days=this_time-start_time
	days=days.days
	return days%21+1
#########################
def SWOT_observation_error():
	"""observation error of WSE depending on the L*W of each pixel
	used sigma*(1/l)*(1/w) l=k*L, w=q*W  Rodrigaz et al 2017:
	According to CaMa k=0.25, q=0.85"""
	# fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	# if os.path.isfile(fname):
	# 	obs_err=np.fromfile(fname,np.float32).reshape(ny,nx)
	# 	return 0
	# else:
	nx,ny,gsize = map_dimension()
	k=1.00 # assume nearest part to the unit catchment
	q=1.00 # used 1.0 -> river width variability is 30%
	ovs_err = 0.10
	rivlen=np.fromfile(CaMa_dir()+"/map/"+mapname()+"/rivlen.bin",np.float32).reshape(ny,nx)
	rivwth=np.fromfile(CaMa_dir()+"/map/"+mapname()+"/rivwth_gwdlr.bin",np.float32).reshape(ny,nx)
	nextx=(np.fromfile(CaMa_dir()+"/map/"+mapname()+"/nextxy.bin",np.int32).reshape(2,ny,nx)[0]!=-9999)*1.0
	rivlen=1.0 #rivlen*1.0e-3 #used as one kilometer
	rivwth=rivwth*1.0e-3
	area=(k*rivlen)*(q*rivwth)
	obs_err=ovs_err*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx
	#obs_err=pm.ovs_err()*(1/(area+1.0e-20))*nextx
	# if water area < 1.0 km2 -> 0.25
	obs_err=obs_err*(area>=1.0)*1.0+0.25*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx*(area<1.0)*1.0
	obs_err=ma.masked_where(area<0.625,obs_err).filled(0.25) # 25cm for < 1km^2 water area
	obs_err=obs_err*((obs_err<=0.25)*1.0) + 0.25*((obs_err>0.25)*1.0)
	obs_err=obs_err*nextx
	obs_err=obs_err.astype(np.float32)
	# obs_err0=obs_err[iy,ix]
	# fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	# obs_err.tofile(fname)
	return obs_err
#########################
def err_rand(obs_err,ix,iy):
	"""make random values to add to true values"""
	# nx,ny,gsize = pm.map_dimension()
	# fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	# obs_err=np.fromfile(fname,np.float32).reshape(ny,nx)
	# obs_err=obs_err*((obs_err<=0.25)*1.0) + 0.25*((obs_err>0.25)*1.0)
	rand = np.random.normal(0.0,obs_err[iy,ix],1)
	return rand
#########################
def dam_loc():
	"Prepare dam location"
	nx,ny,gsize = map_dimension()
	damloc=np.ones([ny,nx],np.float32)
	#-----
	with open(dam_list(),"r") as f:
		lines=f.readlines()
	#
	for line in lines[1::]:
		line    = re.split(" ",line)
		line    = list(filter(None, line))
		#---------------------------------
		ix      = int(line[4]) - 1
		iy      = int(line[5]) - 1
		damloc[iy,ix]=-9999.0
	return damloc
#########################
def write_txt(inputlist):
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    dir0=inputlist[3]
    print ("write text file: ",yyyy,mm,dd)
    target_dt=datetime.date(int(yyyy),int(mm),int(dd))
    txtfile=dir0+"/"+yyyy+mm+dd+".txt"
    if obs_name() == "HydroWeb":
        xlist, ylist, l_wse, m_wse, s_wse, l_sat = HydroWeb_data(yyyy,mm,dd)
    if obs_name() == "SWOT":
        xlist, ylist, l_wse, m_wse, s_wse, l_sat = swot_data(yyyy,mm,dd) 
    #--------------
    pnum=len(xlist)
    # print ('xlist:',pnum, "l_wse:",len(l_wse))
    with open(txtfile,"w") as txtf:
        for point in np.arange(pnum):
            iix=xlist[point]
            iiy=ylist[point]
            wseo=l_wse[point]
            mean_wse=m_wse[point]
            std_wse=s_wse[point]
            sat=l_sat[point]
            line="%04d	%04d	%10.4f	%10.4f	%10.4f	%s\n"%(iix,iiy,wseo,mean_wse,std_wse,sat)
            txtf.write(line)
            print (line)
    return 0
#########################
def prepare_obs(dir0="./"):
	"""
	Prepare observations as textfile
	"""
	# making dir
	mk_dir(dir0)
	#=========================
	syear,smon,sday=starttime()
	eyear,emon,eday=endtime()
	start_dt=datetime.date(syear,smon,sday)
	end_dt=datetime.date(eyear,emon,eday)
	start=0
	last=(end_dt-start_dt).days + 1
	# print (start,last)
	#-------
	inputlist=[]
	for day in np.arange(start,last):
		target_dt=start_dt+datetime.timedelta(days=day)
		yyyy='%04d' % (target_dt.year)
		mm='%02d' % (target_dt.month)
		dd='%02d' % (target_dt.day)
		# print (yyyy,mm,dd) #,obs_dir
		inputlist.append([yyyy,mm,dd,dir0])
	# write text files parallel
	p=Pool(10)
	p.map(write_txt,inputlist)
	p.terminate()
	# map(write_txt,inputlist)
	return 0
####################################
############# parameters ###########
####################################
def starttime():
    return 2002,1,1
####################################
def endtime():
    return 2021,1,1
####################################
def obs_list():
    # return "../dat/HydroWeb_alloc_amz_06min_QC0_simulation.txt"
	# return "../dat/HydroWeb_alloc_amz_06min_2002-2020.txt"
	return "../dat/HydroWeb_alloc_glb_15min.txt"
####################################
def HydroWeb_list():
    # return "../dat/HydroWeb_alloc_amz_06min_QC0_simulation.txt"
	return "../dat/HydroWeb_alloc_amz_06min_2002-2020.txt"
####################################
def obs_name():
    # return "HydroWeb"
	return "SWOT"
####################################
def obs_dir():
    # return "/cluster/data7/menaka/HydroDA/obs/HydroWeb"
    # return "/cluster/data6/menaka/HydroWeb"
    # return "/cluster/data6/menaka/ensemble_org/CaMa_out/E2O003"
	return "/work/a04/julien/CaMa-Flood_v4/out/coupled-model2"
####################################
def dam_list():
	return "../dat/dam_glb_15min.txt"
####################################
def CaMa_dir():
	return "/cluster/data6/menaka/CaMa-Flood_v4"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders
####################################
def mapname():
    # return "amz_06min"
    return "glb_15min"
    # related CaMa-Flood map directory
    # [e.g. : glb_15min, glb_06min, Mkg_06min, etc.]
    # Check 
####################################
def map_dimension():
    fname=CaMa_dir()+"/map/"+mapname()+"/params.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    #-------
    nx     = int(filter(None, re.split(" ",lines[0]))[0])
    ny     = int(filter(None, re.split(" ",lines[1]))[0])
    gsize  = float(filter(None, re.split(" ",lines[3]))[0])
    return nx,ny,gsize
####################################
def out_dir():
	# return "/cluster/data7/menaka/HydroDA/obs/HydroWeb"
	# return "/cluster/data7/menaka/HydroDA/obs/HydroWebAll"
	return "/cluster/data7/menaka/HydroDA/obs/SWOTH08"
####################################
if __name__ == "__main__":
	print ("prepare observations")
	# prepare_obs("/cluster/data7/menaka/HydroDA/obs/HydroWeb")
	# prepare_obs("/cluster/data7/menaka/HydroDA/obs/HydroWeb_glb_15min")
	prepare_obs("/cluster/data7/menaka/HydroDA/obs/SWOTH08") # for SWOTH08
