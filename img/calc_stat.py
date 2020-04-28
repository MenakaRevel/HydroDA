#!/opt/local/bin/python
# -*- coding: utf-8 -*-
""" Calculate statistical maps for river discharge [AI,RMSE, rRMSE, NRMSE, NSE, NSEAI, VE]
 Menaka@IIS 2020/04/22"""
import numpy as np
import datetime
import sys
import os
import string
import calendar
import errno
from numpy import ma
#--
def mk_dir(sdir):
    try:
        os.makedirs(sdir)
    except:
        pass
#---
def slink(src,dst):
  try:
    os.symlink(src,dst)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST:
      os.remove(dst)
      os.symlink(src,dst)
    else:
      raise
#--
def make_yearlist(syear,smon,sday,N=365):
    start_dt=datetime.date(syear,smon,sday)
    fname="year_day.txt"
    f=open(fname,"w")
    for day in np.arange(0,N):
        nw_dt = start_dt + datetime.timedelta(days=day)
        yyyy,mm,dd = nw_dt.year, nw_dt.month, nw_dt.day
        yyyymmdd = "%04d%02d%02d\n"%(yyyy,mm,dd)
        f.write(yyyymmdd)
    f.close()
#--
slink("../gosh/params.py", "params.py")
import params as pm
# define values
N=366 #2004
exp="E2O_wmc_01_anomalyDA"
outdir=pm.DA_dir()+"/out/"+exp
# mkdir stat
mk_dir(outdir+"/assim_out/stat")
# make list of days
#make_yearlist(2004,1,1,N=366)
# compile
#os.system("ifort calc_stat.f90 -o calc_stat -O3 -assume byterecl -heap-arrays  -g -traceback -free -parallel")
# run calc_stat
os.system("./calc_stat "+str(N)+" "+str(pm.ens_mem())+" "+outdir+" "+pm.CaMa_dir())
