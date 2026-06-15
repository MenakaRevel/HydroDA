#! /usr/bin/python

import re
import numpy as np
import math
import os
import sys

from read_patchMS import upstream

sys.path.append("../img")
import read_hydroweb as hweb
#
area_thr = 1.0e09 #m2
slpe_thr = 1.0e20 #m
elev_thr = 1.0e20 #m
dist_thr = 100.0  #km
rmse_thr = 100.0  #m
#-------------------------------------
def slope(ix,iy,nextxy,uparea,elevtn,nxtdst,rivseq):
    nextX=nextxy[0]
    nextY=nextxy[1]
    slp1=0.0
    slp2=0.0
    if rivseq[iy,ix]>1 and nextX[iy,ix]>0:
        uXX, uYY = upstream(ix+1,iy+1,nextX.T,nextY.T,uparea.T)
        uXX = uXX - 1
        uYY = uYY - 1
        slp1=(elevtn[uYY,uXX]-elevtn[iy,ix])
        dXX = nextX[iy,ix] - 1
        dYY = nextY[iy,ix] - 1
        slp2=(elevtn[iy,ix]-elevtn[dYY,dXX])
    elif rivseq[iy,ix]==1:
        slp1=0.0
        dXX = nextX[iy,ix] - 1
        dYY = nextY[iy,ix] - 1
        slp2=(elevtn[iy,ix]-elevtn[dYY,dXX])
    elif nextX[iy,ix]<=0:
        uXX, uYY = upstream(ix+1,iy+1,nextX.T,nextY.T,uparea.T)
        uXX = uXX - 1
        uYY = uYY - 1
        slp1=(elevtn[uYY,uXX]-elevtn[iy,ix])
        slp2=0.0
    return slp1,slp2
#===================
syear=2002
<<<<<<< HEAD
eyear=2014
#===================
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
CaMa_dir="/cluster/data7/menaka/CaMa-Flood_v407"
map="glb_15min"
=======
eyear=2020
#===================
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# map="glb_15min"
>>>>>>> dev_virtual
# map="glb_06min"
# map="amz_06min"
fname=CaMa_dir+"/map/"+map+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
# gsize=0.25
# nx=1440
# ny=720
nextxy = CaMa_dir+"/map/"+map+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+map+"/rivwth.bin"
rivhgt = CaMa_dir+"/map/"+map+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+map+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+map+"/elevtn.bin"
uparea = CaMa_dir+"/map/"+map+"/uparea.bin"
lonlat = CaMa_dir+"/map/"+map+"/lonlat.bin"
nxtdst = CaMa_dir+"/map/"+map+"/nxtdst.bin"
rivseq = CaMa_dir+"/map/"+map+"/rivseq.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#---
clas =  CaMa_dir+"/map/"+map+"/class.bin"
# clas =  np.fromfile(clas,np.int32).reshape(ny,nx)
#----
dx=5
dy=5
#--
#fname="HydroWebStation_MERITv07.csv"
# fname="./data/HydroWeb_VS"
# f=open(fname,"r")
# lines=f.readlines()
# f.close()
# write
#---------------------
# rivername0="AMAZONAS" #"CONGO" #
rivername0="ALL" #"CONGO" #
stream0=["AMAZONAS","SOLIMOES"] # "CONGO" #
#---------------------
# fname="HydroWeb_alloc_"+map+".txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+map+"_20210826.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+map+"_20210826.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+map+"_20210920.txt"
fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+map+"_20220725.txt"
with open(fname,"r") as f:
	lines=f.readlines()
#===============================================
# unreal observations
# fname="/cluster/data6/menaka/Altimetry/out/unreal_obs_20210824.txt"
# fname="/cluster/data6/menaka/Altimetry/out/unreal_obs_20210826.txt"
fname="/cluster/data6/menaka/AltiMaP/out/unreal_obs_20220802.txt"
with open(fname,"r") as f_unreal:
	unreal=f_unreal.readlines()
#----
unreal_stations=[]
for item in unreal:
    item    = list(filter(None, re.split(" ",item)))
    station = item[0]
    unreal_stations.append(station)
#===============================================
# higher RMSE locations
# fname="/cluster/data6/menaka/AltiMaP/out/rmse_VIC_BC_glb_06min_20211122.txt"
fname="/cluster/data6/menaka/AltiMaP/out/rmse_VIC_BC_glb_15min_20221020.txt"
with open(fname,"r") as f_rmse:
	rmse=f_rmse.readlines()
#----
rmse_stations=[]
for item in rmse[1::]:
    item    = list(filter(None, re.split(" ",item)))
    station = item[0]
    rmse0   = float(item[1])
    if rmse0 > rmse_thr:
        rmse_stations.append(station)
#===
# writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+".txt"
# writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_amz.txt"
# writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QC.txt"
# writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QC1.txt"
# writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QC0.txt"
# writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QCrmse.txt"
<<<<<<< HEAD
writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+".txt"
=======
writef="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_2002-2020.txt"
>>>>>>> dev_virtual
with open(writef, "w") as wf:
    header = "%13s%62s%8s%8s%8s%8s%12s%12s%12s%12s%17s\n"%("ID","station","lon","lat","ix","iy","elevation","ele_diff","EGM08","EGM96","satellite")
    wf.write(header)
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        #print line
        num     = int(line[0])
        station = line[1]
        line2   = re.split("_",station)
        riv     = line2[1]
        stream  = line2[2]
        lon     = float(line[3])
        lat     = float(line[4])
        ix      = int(line[5])
        iy      = int(line[6])
        ele     = elevtn[iy-1,ix-1] #float(line[7])
        ele_dif = float(line[7])
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        sat     = line[10].strip()
        dist    = float(line[11])
        flag    = int(line[12])
        kx      = int(line[13])
        ky      = int(line[14])
        #======================
        if rivername0 != "ALL":
            if riv != rivername0:
                continue
        # if stream not in stream0:
        #     continue
        # print (num, ele_dif)

        ################
        # data available years
        ################
        org=hweb.HydroWeb_continous_WSE(station,syear=syear,eyear=eyear)
        if np.sum((org!=-9999.0)*1.0) == 0.0:
            print ("no data available for ",syear,"-",eyear)
            continue

        ################
        # condtion for elevation differnce
        ################
        if abs(ele_dif) > elev_thr:
            print ("elevation differnce is too large: ", ele_dif)
            continue

        ################
        # condtion for mainstream
        ################
        if uparea[iy-1,ix-1] < area_thr:
            print ("smaller river: ",station, uparea[iy-1,ix-1])
            continue
        
        ################
        # condition for slope
        ################
        slp1, slp2 =slope(ix-1,iy-1,nextxy,uparea,elevtn,nxtdst,rivseq)
        if slp1 > slpe_thr:
            print ("high slope:",station, slp1)
            continue

        ################
        # condition for distance to mouth
        ################
        if dist > dist_thr:
            print ("large dist to mouth:",station, dist)
            continue

        ################
        # condition for unreal observations
        ################
        if station in unreal_stations:
            print ("unreal observations: ", station)
            continue

        ################
        # condition for higher rmse observations
        ################
        if station in rmse_stations:
            print ("hihger RMSE virtual station: ( >",rmse_thr,"m)", station)
            continue

        linew="%013d%62s%8.2f%8.2f%8d%8d%12.2f%12.2f%12.2f%12.2f%17s\n"%(num,station,lon,lat,ix,iy,ele,ele_dif,EGM08,EGM96,sat)
        # print (linew)
        wf.write(linew)