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
import json

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
####################################
# SWOT
#########################
def swot_data(yyyy,mm,dd):
	# prepare synthetic observations using
	# pre-simulated data
	# river width thershold
	obs=obs_SWOT(yyyy,mm,dd)
	# rivwdth_thr=50.0 #m
	nx,ny,gsize = map_dimension()
	# ny_swot = min(ny,640)
	# day=SWOT_day(yyyy,mm,dd)
	# SWOTDD="%02d"%(day)
	# fname="../sat/mesh_day"+SWOTDD+".bin" # for glb_15min
	# mesh_in=np.fromfile(fname,np.float32).reshape([ny_swot,nx])
	# mesh=(mesh_in>=10)*(mesh_in<=60)
	# meshP=mesh-1000*(mesh<0.1)
	# SWOTmesh=np.zeros([ny,nx],np.float32)
	# SWOTmesh[40:680,:]=meshP
	# fname=CaMa_dir()+"/map/"+mapname()+"/rivwth_gwdlr.bin"
	# rivwth=np.fromfile(fname,np.float32).reshape(ny,nx)
	# damloc=dam_loc()
	# obs=(SWOTmesh>=1.0)*(rivwth>=rivwdth_thr)*(damloc>0.0)*1.0
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
	nt = 366 if calendar.isleap(year) else 365
	# if calendar.isleap(year):
	# 	nt=366
	# else:
	# 	nt=365
	orgfile=np.fromfile(fname,np.float32).reshape([nt,ny,nx])
	#===================================
	date = datetime.date(year, mon, day)
	day_of_year = date.timetuple().tm_yday
	seed=int(day_of_year)
	obs_err=SWOT_observation_error(seed)
	#-----------------------
	start_dt=datetime.date(year,1,1)
	target_dt=datetime.date(year,mon,day)
	it=(target_dt-start_dt).days
	#-----------------------
	# calculate mean and std of wse
	syear,smon,sday=starttime()
	eyear,emon,eday=endtime()
	SWOTmean,SWOTstd=calc_stat_SWOT(syear,eyear)
	for ix in np.arange(nx):
		for iy in np.arange(ny):
			# print (obs[iy,ix])
			if obs[iy,ix] == 1.0:
				random.seed(seed)
				wse=orgfile[it,iy,ix] + err_rand(obs_err,ix,iy,seed)
				# print (ix,iy,wse[0])
				l_wse.append(wse[0])
				xlist.append(ix+1)
				ylist.append(iy+1)
				m_wse.append(SWOTmean[iy,ix])
				s_wse.append(SWOTstd[iy,ix])
				l_sat.append("SWOT")
	return xlist, ylist, l_wse, m_wse, s_wse, l_sat
####################################
def calc_stat_SWOT(syear,eyear,obs_dir):
	"""calculate statistics of SWOT data"""
	nx,ny,gsize = map_dimension()
	ny_swot = min(ny,640)
	# obs=np.zeros([get_days(syear, eyear),ny,nx],np.float32)
	SWOTOBS=np.zeros([get_days(syear, eyear),ny,nx],np.float32)
	indays=0
	for year in np.arange(syear,eyear+1):
		nt = 366 if calendar.isleap(year) else 365
		#-----------------------
		start_dt=datetime.date(year,1,1)
		target_dt=datetime.date(year,12,31)
		it=(target_dt-start_dt).days
		#-----------------------
		odir=obs_dir
		fname=odir+"/sfcelv"+str(year)+".bin"
		SWOTOBS[indays:indays+get_days(syear,year),:,:]=np.fromfile(fname,np.float32).reshape([nt,ny,nx])
		print (indays,get_days(syear,year))
		indays=indays+get_days(syear,year)
		# obs_err=SWOT_observation_error()
		#-----------------------
		for day in np.arange(nt):
			target_dt=start_dt+datetime.timedelta(days=day)
			yyyy="%04d"%(target_dt.year)
			mm="%02d"%(target_dt.month)
			dd="%02d"%(target_dt.day)
			#-----------------------
			random.seed(day)
			obs_err=SWOT_observation_error(day)
			# print (day, obs_SWOT(yyyy,mm,dd)[100,100]),err_rand_array(obs_err,day)[100,100]
			# obs[day,:,:]=obs_SWOT(yyyy,mm,dd)
			SWOTOBS[day,:,:]=SWOTOBS[day,:,:]+err_rand_array(obs_err,day)
			obs=obs_SWOT(yyyy,mm,dd)
			SWOTOBS[day,:,:]=ma.masked_where(obs!=1.0,SWOTOBS[day,:,:]).filled(-9999.0)
			# print (day, SWOTOBS[day,100,100])
	#-----------------------
	# SWOTOBS=ma.masked_where(obs!=1.0,orgfile).filled(-9999.0)
	# calculate statistics
	#-----------------------
	# mean
	mean=np.mean(ma.masked_equal(SWOTOBS,-9999.0),axis=0)
	std=np.std(ma.masked_equal(SWOTOBS,-9999.0),axis=0)
	#-----------------------
	return SWOTOBS,mean,std
####################################
def obs_SWOT(yyyy,mm,dd):
	"""get SWOT observation days"""
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
	return obs
####################################
def SWOT_day(yyyy,mm,dd):
	st_year,st_month,st_date=starttime()
	start_time=datetime.date(st_year,st_month,st_date)
	this_time=datetime.date(int(yyyy),int(mm),int(dd))
	days=this_time-start_time
	days=days.days
	return days%21+1
#########################
def SWOT_observation_error(seed):
	"""observation error of WSE depending on the L*W of each pixel
	used sigma*(1/l)*(1/w) l=k*L, w=q*W  Rodrigaz et al 2017:
	According to CaMa k=0.25, q=0.85"""
	# fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	# if os.path.isfile(fname):
	# 	obs_err=np.fromfile(fname,np.float32).reshape(ny,nx)
	# 	return 0
	# else:
	random.seed(seed)
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
def err_rand(obs_err,ix,iy,seed):
	"""make random values to add to true values"""
	# nx,ny,gsize = pm.map_dimension()
	# fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
	# obs_err=np.fromfile(fname,np.float32).reshape(ny,nx)
	# obs_err=obs_err*((obs_err<=0.25)*1.0) + 0.25*((obs_err>0.25)*1.0)
	random.seed(seed)
	rand = np.random.normal(0.0,obs_err[iy,ix],1)
	return rand
#########################
def err_rand_array(obs_err,seed):
    """make random values to add to true values"""
    # nx,ny,gsize = pm.map_dimension()
    # fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/obs/obs_err.bin"
    # obs_err=np.fromfile(fname,np.float32).reshape(ny,nx)
    # obs_err=obs_err*((obs_err<=0.25)*1.0) + 0.25*((obs_err>0.25)*1.0)
    ny,nx=obs_err.shape
    # obs_err1=np.zeros((ny,nx),np.float32)
    # obs_err1 = [[(random.seed(seed+iy*nx+ix), np.random.normal(0.0, obs_err[iy, ix], 1))[1] for ix in range(nx)] for iy in range(ny)]
    # obs_err1 = np.array(obs_err1).reshape(ny,nx)

    # Generate random seeds for each element in the array
    # np.random.seed(seed + np.arange(ny)[:, np.newaxis] * nx + np.arange(nx))
    # Generate random values using vectorized operations
    obs_err1 = np.random.normal(0.0, obs_err, (ny, nx))
	# for iy in range(ny):
	# 	for ix in range(nx):
	# 		seed=seed+iy*nx+ix
	# 		random.seed(seed)
	# 		rand = np.random.normal(0.0,obs_err[iy,ix],1)
	# 		obs_err1[iy,ix]=rand
    return obs_err1
#########################
def get_days(syear, eyear):
    start_date = datetime.date(syear, 1, 1)
    end_date = datetime.date(eyear, 12, 31)
    delta = (end_date - start_date).days + 1
    return delta
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
    xlist, ylist, l_wse, m_wse, s_wse, l_sat = inputlist
    print ("write text file: ",yyyy,mm,dd)
    target_dt=datetime.date(int(yyyy),int(mm),int(dd))
    txtfile=dir0+"/"+yyyy+mm+dd+".txt"
    # if obs_name() == "HydroWeb":
    # 	xlist, ylist, l_wse, m_wse, s_wse, l_sat = HydroWeb_data(yyyy,mm,dd)
    # if obs_name() == "SWOT":
    # 	xlist, ylist, l_wse, m_wse, s_wse, l_sat = swot_data(yyyy,mm,dd) 
    # if obs_name() == "CGLS":
    # 	xlist, ylist, l_wse, m_wse, s_wse, l_sat = CGLS_data(yyyy,mm,dd)
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
def prepare_obs(inputlist):
	"""
	Prepare observations as textfile
	"""
	indir=inputlist[0]
	outdir=inputlist[1]
	# making dir
	mk_dir(outdir)
	#=========================
	nx,ny,gsize = map_dimension()
	syear,smon,sday=starttime()
	eyear,emon,eday=endtime()
	start_dt=datetime.date(syear,smon,sday)
	end_dt=datetime.date(eyear,emon,eday)
	start=0
	last=(end_dt-start_dt).days + 1
	#--
	print ("Calculation of statistics of observations")
	SWOTobs,SWOTmean,SWOTstd=calc_stat_SWOT(syear,eyear,indir)
	#--------------
	for day in np.arange(start,last):
		target_dt=start_dt+datetime.timedelta(days=day)
		yyyy='%04d' % (target_dt.year)
		mm='%02d' % (target_dt.month)
		dd='%02d' % (target_dt.day)
		txtfile=outdir+"/"+yyyy+mm+dd+".txt"
		print ("write text file: ",yyyy,mm,dd)
		with open(txtfile,"w") as txtf:
			for ix in np.arange(nx):
				for iy in np.arange(ny):
					# print (ix,iy,SWOTobs[day,iy,ix],SWOTmean[iy,ix],SWOTstd[iy,ix])
					if SWOTobs[day,iy,ix] != -9999.0:
						line="%04d	%04d	%10.4f	%10.4f	%10.4f	%s\n"%(ix,iy,SWOTobs[day,iy,ix],SWOTmean[iy,ix],SWOTstd[iy,ix],"SWOT")
						txtf.write(line)
						print (line)
	return 0
####################################
############# parameters ###########
####################################
def starttime():
    return 2001,1,1
####################################
def endtime():
    return 2001,12,31
####################################
def obs_list(): ## only for real observations
    # return "../dat/HydroWeb_alloc_amz_06min_QC0_simulation.txt"
	# return "../dat/HydroWeb_alloc_amz_06min_2002-2020.txt"
	# return "../dat/HydroWeb_alloc_glb_15min.txt"
	# return "../dat/HydroWeb_alloc_conus_06min_org.txt"
	# return "../dat/HydroWeb_alloc_conus_06min_DIR.txt"
	# return "../dat/CGLS_alloc_conus_06min_DIR.txt"
	return "../dat/CGLS_alloc_conus_06min_org.txt"
####################################
def HydroWeb_list(): ### not used 
    # return "../dat/HydroWeb_alloc_amz_06min_QC0_simulation.txt"
	# return "../dat/HydroWeb_alloc_amz_06min_2002-2020.txt"
	# return "../dat/HydroWeb_alloc_conus_06min_org.txt"
	# return "../dat/HydroWeb_alloc_conus_06min_DIR.txt"
	# return "../dat/CGLS_alloc_conus_06min_DIR.txt"
	return "../dat/CGLS_alloc_conus_06min_org.txt"
####################################
def obs_name():
    # return "HydroWeb"
	return "SWOT"
	# return "CGLS"
####################################
def obs_dir():
    # return "/cluster/data7/menaka/HydroDA/obs/HydroWeb"
    # return "/cluster/data6/menaka/HydroWeb"
	# return "/work/a06/menaka/CGLS"
    # return "/cluster/data6/menaka/ensemble_org/CaMa_out/E2O003"
	# return "/work/a04/julien/CaMa-Flood_v4/out/coupled-model2"
	# return "/cluster/data6/menaka/CaMa-H08/out/obs_org"
	# return "/cluster/data6/menaka/CaMa-H08/out/obs_rivhgt"
	# return "/cluster/data6/menaka/CaMa-H08/out/obs_rivwth"
	# return "/cluster/data6/menaka/CaMa-H08/out/obs_rivman"
	# return "/cluster/data6/menaka/CaMa-H08/out/obs_fldhgt"
	return "/cluster/data6/menaka/CaMa-H08/out/obs_corr_all_001"
####################################
def dam_list():
	return "../dat/dam_glb_15min.txt"
	# return "../dat/dam_conus_06min.txt"
####################################
def CaMa_dir():
	return "/cluster/data6/menaka/CaMa-Flood_v4"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders
####################################
def mapname():
    # return "amz_06min"
    return "glb_15min"
	# return "conus_06min"
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
def out_dir(): ### not used --> give as direct input @L536
	# return "/cluster/data7/menaka/HydroDA/obs/HydroWeb"
	# return "/cluster/data7/menaka/HydroDA/obs/HydroWebAll"
	# return "/cluster/data7/menaka/HydroDA/obs/SWOTH08"
	# return "/cluster/data7/menaka/HydroDA/obs/SWOT_CaMaH08_org"
	# return "/cluster/data7/menaka/HydroDA/obs/SWOT_CaMaH08_rivhgt"
	# return "/cluster/data7/menaka/HydroDA/obs/SWOT_CaMaH08_fldhgt"
	return "/cluster/data7/menaka/HydroDA/obs/SWOT_CaMaH08_all"
####################################
#########################
if __name__=="__main__":
	if len(sys.argv) > 1:
		ncpus=int(sys.argv[1])
	else:
		ncpus=20
	indir0="/cluster/data6/menaka/CaMa-H08/out"
	outdir0="/cluster/data7/menaka/HydroDA/obs"
	inputlist=[]
	#--------------
	for exp in range(1,20+1):
		indir=indir0+"/obs_corr_all_%03d"%(exp)
		outdir=outdir0+"/SWOT_CaMaH08_all_%03d"%(exp)
		inputlist.append([indir,outdir])
	#--------------
	#prepare SWOT observations parallel
	p=Pool(ncpus)
	p.map(prepare_obs,inputlist)
	p.terminate()
    
    # indir=indir0+"/obs_corr_all_001"
    # outdir=outdir0+"/SWOT_CaMaH08_all_001"
    # prepare_obs([indir,outdir])
	# print ("Calculation of statistics")
	# SWOTobs,SWOTmean,SWOTstd=calc_stat_SWOT(2001,2001)