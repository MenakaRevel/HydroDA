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
slink("../gosh/params.py", "params.py")
import params as pm
# define values
N=366 #2004
exp="E2O_womc"
outdir=pm.HydroDA()+"/out/"+exp
# mkdir stat
mk_dir(outdir+"/assim_out/stat")
# compile
os.system("ifort calc_stat.f90 -o calc_stat -O3 -assume byterecl -heap-arrays  -g -traceback -free -parallel")
# run calc_stat
os.system("./calc_stat "+str(N)+" "+str(pm.ens_mem())+" "+out_dir+" "+pm.CaMa_dir())
