#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
to read and write SWOT observations
'''
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
import pandas as pd
import geopandas
import requests
from io import StringIO
#################################################################
#################################################################
#=========================================================================================
# Function 
#=========================================================================================
###########################
def mk_dir(sdir):
    try:
        os.makedirs(sdir)
    except:
        pass
###########################
def obs_list():
    # return "../dat/HydroWeb_alloc_amz_06min_QC0_simulation.txt"
	# return "../dat/HydroWeb_alloc_amz_06min_2002-2020.txt"
	# return "../dat/HydroWeb_alloc_glb_15min.txt"
	# return "../dat/HydroWeb_alloc_conus_06min_org.txt"
	# return "../dat/HydroWeb_alloc_conus_06min_DIR.txt"
	# return "../dat/CGLS_alloc_conus_06min_DIR.txt"
	# return "../dat/CGLS_alloc_conus_06min_org.txt"
    return "../dat/SWOT_alloc_Mackenzie_06min.txt"
###########################
def query_hydrocron(station, ix=-9999, iy=-9999, EGM08=0.0, EGM96=0.0, start_time="2024-01-01T00:00:00Z", end_time="2024-11-01T00:00:00Z"):
    """Query Hydrocron for reach-level time series data.

    Parameters
    ----------
    station: str - String SWORD reach identifier
    start_time: str - String time to start query format: "2024-01-01T00:00:00Z"
    end_time: str - String time to end query

    Returns
    -------
    pandas.DataFrame that contains query results
    """
    query_url= "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"
    fields = "node_id,time_str,wse,wse_u,wse_r_u,node_q"
    params = {
        "feature": 'Node',
        "feature_id": station,
        "output": "csv",
        "start_time": start_time,
        "end_time": end_time,
        "fields": fields
    }
    results = requests.get(query_url, params=params)
    if "results" in results.json().keys():
        results_csv = results.json()["results"]["csv"]
        df = pd.read_csv(StringIO(results_csv))
        # Remove fill values for missing observations
        df = df.loc[(df["wse"] != -999999999999.0)].reset_index(drop=True)
        df = df.loc[(df['node_q']<=1) & (df['wse_r_u']<1.0), :]
        # Convert time_str to datetime format
        df.time_str = pd.to_datetime(df.time_str)
        # add x and y
        df['x'] = ix
        df['y'] = iy
        df['EGM08'] = EGM08
        df['EGM96'] = EGM96
        df = df.loc[:,['node_id','time_str','x','y','wse','EGM08','EGM96']]
    else:
        df = pd.DataFrame({
            "node_id": np.int64(station),
            "time_str": datetime.datetime(1900, 1, 1).strftime("%Y-%m-%dT%H:%M:%S"),
            "x":-9999,
            "y":-9999,
            "wse": -999999999999.0,
            "EGM08": -999999999999.0,
            "EGM96": -999999999999.0,
            }, index=[0])

    return df
####################################
def starttime():
    return 2024,1,1
####################################
def endtime():
    return 2024,11,1
# #########################
# def write_txt(inputlist):
# 	yyyy=inputlist[0]
# 	mm=inputlist[1]
# 	dd=inputlist[2]
# 	dir0=inputlist[3]
# 	print ("write text file: ",yyyy,mm,dd)
# 	target_dt=datetime.date(int(yyyy),int(mm),int(dd))
# 	txtfile=dir0+"/"+yyyy+mm+dd+".txt"
#     lstan, xcods, ycods, leledif, lEGM08, lEGM96, satellite = get_HydroWeb() #get_CGLS()
# 	pnum=len(lstan)
#     lstan[point],lEGM08[point],lEGM96[point]

#     # df_w=ddf[ddf['node_id']==]

#     #--------------
# 	pnum=len(xlist)
# 	# print ('xlist:',pnum, "l_wse:",len(l_wse))
# 	with open(txtfile,"w") as txtf:
# 		for point in np.arange(pnum):
# 			iix=xlist[point]
# 			iiy=ylist[point]
# 			wseo=l_wse[point]
# 			mean_wse=m_wse[point]
# 			std_wse=s_wse[point]
# 			sat=l_sat[point]
# 			line="%04d	%04d	%10.4f	%10.4f	%10.4f	%s\n"%(iix,iiy,wseo,mean_wse,std_wse,sat)
# 			txtf.write(line)
# 			print (line)
# 	return 0
# #########################
# def prepare_obs(dir0="./"):
# 	"""
# 	Prepare observations as textfile
# 	"""
# 	# making dir
# 	mk_dir(dir0)
# 	#=========================
# 	syear,smon,sday=starttime()
# 	eyear,emon,eday=endtime()
# 	start_dt=datetime.date(syear,smon,sday)
# 	end_dt=datetime.date(eyear,emon,eday)
# 	start=0
# 	last=(end_dt-start_dt).days + 1
# 	# print (start,last)
# 	#-------
# 	inputlist=[]
# 	for day in np.arange(start,last):
# 		target_dt=start_dt+datetime.timedelta(days=day)
# 		yyyy='%04d' % (target_dt.year)
# 		mm='%02d' % (target_dt.month)
# 		dd='%02d' % (target_dt.day)
# 		# print (yyyy,mm,dd) #,obs_dir
# 		inputlist.append([yyyy,mm,dd,dir0])
# 	# write text files parallel
# 	p=Pool(20)
# 	p.map(write_txt,inputlist)
# 	p.terminate()
# 	# map(write_txt,inputlist)
# 	return 0
#########################
results = []
fname=obs_list()
# fname=HydroWeb_list()
#=========================
with open(fname,"r") as f:
    lines=f.readlines()
for line in lines[1::]:
    line    = re.split(" ",line)
    line    = list(filter(None, line))
    #print line
    num     = line[0]
    station = line[1]
    # riv     = re.split("_",station)[1]
    # lon     = float(line[2])
    # lat     = float(line[3])
    ix      = int(line[4])
    iy      = int(line[5])
    eledif  = float(line[7])
    EGM08   = float(line[8])
    EGM96   = float(line[9])
    sat     = line[10].split()[0]
    print (station)
    results.append(query_hydrocron(station, ix=ix, iy=iy, EGM08=EGM08, EGM96=EGM96))
#======================================
# Load DataFrame results into dataframe
ddf = pd.concat(results)
ddf.head(n=20)
# Remove fill values for missing observations
ddf = ddf.loc[(ddf["wse"] != -999999999999.0)]

# Convert time_str to datetime format
ddf.time_str = pd.to_datetime(ddf.time_str)

print (ddf)
print (ddf['time_str'].dt.date.unique())
#===
# making dir
dir0='/cluster/data7/menaka/HydroDA/obs/SWOT_Mackenzie_06min'
mk_dir(dir0)
#=========================
syear,smon,sday=starttime()
eyear,emon,eday=endtime()
start_dt=datetime.date(syear,smon,sday)
end_dt=datetime.date(eyear,emon,eday)
start=0
last=(end_dt-start_dt).days + 1
for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=int(day))
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    #===========================
    print ("write text file: ",yyyy,mm,dd)
    target_dt=datetime.date(int(yyyy),int(mm),int(dd))
    txtfile=dir0+"/"+yyyy+mm+dd+".txt"
    #===========================
    with open(txtfile,"w") as txtf:
        for station_id in ddf['node_id'].dropna().unique():
            df_wse   = ddf[ddf['node_id']==station_id]
            iix      = df_wse['x'].values[0]
            iiy      = df_wse['y'].values[0]
            # print (df_wse)
            if not (df_wse['time_str'].dt.date==target_dt).any():
                print ("no obs : ", target_dt.strftime('%Y-%m-%d'), station_id)
                continue
            wseo     = df_wse.loc[df_wse['time_str'].dt.date==target_dt,'wse'].values[0]
            if wseo < -9999.0:
                continue
            mean_wse = df_wse['wse'].mean()
            std_wse  = df_wse['wse'].std()
            sat      = 'SWOT'
            line="%04d	%04d	%10.4f	%10.4f	%10.4f	%s\n"%(iix,iiy,wseo,mean_wse,std_wse,sat)
            txtf.write(line)
            print (line)
