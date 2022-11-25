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
	fname=pm.obs_list()
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
	nx,ny,gsize = pm.map_dimension()
	ny_swot = min(ny,640)
	day=SWOT_day(yyyy,mm,dd)
	SWOTDD="%02d"%(day)
	fname="../../sat/mesh_day"+SWOTDD+".bin" # for glb_15min
	mesh_in=np.fromfile(fname,np.float32).reshape([ny_swot,nx])
	mesh=(mesh_in>=10)*(mesh_in<=60)
	meshP=mesh-1000*(mesh<0.1)
	SWOTmesh=np.zeros([ny,nx],np.float32)
	SWOTmesh[40:680,:]=meshP
	fname=pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivwth_gwdlr.bin"
	rivwth=np.fromfile(fname,np.float32).reshape(ny,nx)
	obs=(SWOTmesh>=1.0)*(rivwth>=rivwdth_thr)*1.0
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
	odir=pm.obs_dir()
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
				wse=orgfile[it,iy,ix] + err_rand(ix,iy)
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
	st_year,st_month,st_date=pm.starttime()
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
	fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	if os.path.isfile(fname):
		return 0
	else:
		nx,ny,gsize = pm.map_dimension()
		k=1.00 # assume nearest part to the unit catchment
		q=1.00 # used 1.0 -> river width variability is 30%
		ovs_err = 0.10
		rivlen=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivlen.bin",np.float32).reshape(ny,nx)
		rivwth=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivwth_gwdlr.bin",np.float32).reshape(ny,nx)
		nextx=(np.fromfile(pm.CaMa_dir()+"/map/glb_15min/nextxy.bin",np.int32).reshape(2,ny,nx)[0]!=-9999)*1.0
		rivlen=1.0 #rivlen*1.0e-3 #used as one kilometer
		rivwth=rivwth*1.0e-3
		area=(k*rivlen)*(q*rivwth)
		obs_err=ovs_err*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx
		#obs_err=pm.ovs_err()*(1/(area+1.0e-20))*nextx
		# if water area < 1.0 km2 -> 0.25
		obs_err=obs_err*(area>=1.0)*1.0+0.25*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx*(area<1.0)*1.0
		obs_err=ma.masked_where(area<0.625,obs_err).filled(0.25) # 25cm for < 1km^2 water area
		obs_err=obs_err*nextx
		obs_err=obs_err.astype(np.float32)
		# obs_err0=obs_err[iy,ix]
		fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
		obs_err.tofile(fname)
		return 0
#########################
def err_rand(ix,iy):
	"""make random values to add to true values"""
	nx,ny,gsize = pm.map_dimension()
	fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	obs_err=np.fromfile(fname,np.float32).reshape(ny,nx)
	obs_err=obs_err*((obs_err<=0.25)*1.0) + 0.25*((obs_err>0.25)*1.0)
	rand = np.random.normal(0.0,obs_err[iy,ix],1)
	return rand
#########################
def write_txt(inputlist):
	yyyy=inputlist[0]
	mm=inputlist[1]
	dd=inputlist[2]
	print ("write text file: ",yyyy,mm,dd)
	# obs_dir="/cluster/data6/menaka/HydroWeb"
	# if pm.obs_name() == "HydroWeb":
	# 	obs_dir="/cluster/data6/menaka/HydroWeb"
	# HydroWeb_dir=inputlist[3]
	target_dt=datetime.date(int(yyyy),int(mm),int(dd))
	txtfile=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/"+yyyy+mm+dd+".txt"
	# print txtfile
	# pnum=len(lname)
	# print (pnum)
	# print pnum
	# == read relevant observation data ==
	# print (yyyy,mm,dd,pm.obs_name())
	if pm.obs_name() == "HydroWeb":
		xlist, ylist, l_wse, m_wse, s_wse, l_sat = HydroWeb_data(yyyy,mm,dd)
	if pm.obs_name() == "SWOT":
		xlist, ylist, l_wse, m_wse, s_wse, l_sat = swot_data(yyyy,mm,dd) 
	#--------------
	pnum=len(xlist)
	# print ('xlist:',pnum, "l_wse:",len(l_wse))
	with open(txtfile,"w") as txtf:
		# for point in np.arange(pnum):
		# 	# == read relevant observation data ==
		# 	if pm.obs_name() == "HydroWeb":
		# 		# == for HydroWeb data ==
		# 		wseo, mean_wse, std_wse = read_HydroWeb(yyyy,mm,dd,lname[point],lEGM08[point],lEGM96[point],leledif[point])
		# 	else:
		# 		wseo, mean_wse, std_wse = read_HydroWeb(yyyy,mm,dd,lname[point],lEGM08[point],lEGM96[point],leledif[point])
		# 	print (point, lname[point], wseo)
		# 	if wseo == -9999.0:
		# 		continue
		# 	iix=xlist[point]
		# 	iiy=ylist[point]
			# mean_wse=np.mean(np.array(lwse))
			# std_wse=np.std(np.array(lwse))
			# sat=satellite[point]
		# pnum=len(xlist)
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
def prepare_obs_old():
	"""
	Prepare observations as textfile
	"""
	# global lname, xlist, ylist, leledif, lEGM08, lEGM96, satellite
	# if pm.obs_name() == "HydroWeb":
	# 	lname, xlist, ylist, leledif, lEGM08, lEGM96, satellite = get_HydroWeb()
	# 	print ("lname: ",len(lname))
	# else:
	# 	lname, xlist, ylist, leledif, lEGM08, lEGM96, satellite = get_HydroWeb()
	# print (len(lname))
	syear,smon,sday=pm.starttime()
	eyear,emon,eday=pm.endtime()
	start_dt=datetime.date(syear,smon,sday)
	end_dt=datetime.date(eyear,emon,eday)
	start=0
	last=(end_dt-start_dt).days
	# print (start,last)
	#-------
	inputlist=[]
	for day in np.arange(start,last):
		target_dt=start_dt+datetime.timedelta(days=day)
		yyyy='%04d' % (target_dt.year)
		mm='%02d' % (target_dt.month)
		dd='%02d' % (target_dt.day)
		# print (yyyy,mm,dd) #,obs_dir
		inputlist.append([yyyy,mm,dd])
	# write text files parallel
	p=Pool(pm.cpu_nums()*pm.para_nums())
	p.map(write_txt,inputlist)
	p.terminate()
	# map(write_txt,inputlist)
	return 0
####################################
def prepare_obs():
	"""
	link observation files
	"""
	if os.path.islink("./assim_out/obs"):
		os.system("rm -rf ./assim_out/obs")
	os.system("ln -sf "+pm.obs_dir()+" ./assim_out/obs")
	return 0
####################################
if __name__ == "__main__":
	print ("prepare observations")
	prepare_obs()