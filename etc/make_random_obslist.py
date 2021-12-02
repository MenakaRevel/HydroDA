#! /usr/bin/python

import re
import numpy as np
import math
import random
#==============
simprt=0.80
#==============
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
#map="glb_15min"
# map="glb_06min"
map="amz_06min"
#==============
readf="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QC0.txt"
# readf="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QCrmse.txt"
#===============================
with open(readf, "r") as rf:
    lines = rf.readlines()
station_list=[]
for line in lines[1::]:
    line    = re.split(" ",line)
    line    = list(filter(None, line))
    #print line
    num     = int(line[0])
    station = line[1]
    station_list.append(station)
#===============================
pnum=len(station_list)
simnum=int(pnum*simprt)
valnum=pnum-simnum
simstation_list=random.sample(station_list,simnum)
#===============================================
# list of virtual stations for data assimilation
simf="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QC0_simulation.txt"
with open(simf, "w") as sf:
    header = "%13s%62s%8s%8s%8s%8s%12s%12s%12s%12s%17s\n"%("ID","station","lon","lat","ix","iy","elevation","ele_diff","EGM08","EGM96","satellite")
    sf.write(header)
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        num     = int(line[0])
        station = line[1]
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        ele     = float(line[6])
        ele_dif = float(line[7])
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        sat     = line[10].strip()
        #======================
        if station in simstation_list:
            linew="%013d%62s%8.2f%8.2f%8d%8d%12.2f%12.2f%12.2f%12.2f%17s\n"%(num,station,lon,lat,ix,iy,ele,ele_dif,EGM08,EGM96,sat)
            # print (linew)
            sf.write(linew)
#========================================
# list of virtual stations for validation
valf="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+map+"_QC0_validation.txt"
with open(valf, "w") as vf:
    header = "%13s%62s%8s%8s%8s%8s%12s%12s%12s%12s%17s\n"%("ID","station","lon","lat","ix","iy","elevation","ele_diff","EGM08","EGM96","satellite")
    vf.write(header)
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        #print line
        num     = int(line[0])
        station = line[1]
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        ele     = float(line[6])
        ele_dif = float(line[7])
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        sat     = line[10].strip()
        #======================
        if station not in simstation_list:
            linew="%013d%62s%8.2f%8.2f%8d%8d%12.2f%12.2f%12.2f%12.2f%17s\n"%(num,station,lon,lat,ix,iy,ele,ele_dif,EGM08,EGM96,sat)
            # print (linew)
            vf.write(linew)