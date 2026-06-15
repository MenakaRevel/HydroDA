#! /usr/bin/python

import re
import numpy as np
import math
import os
import sys

sys.path.append("../img")
import read_grdc as grdc
#
thr_yr=0 #years
#===================
syear=2009
eyear=2014
#===================
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
#map="glb_15min"
# map="glb_06min"
map="amz_06min"

fname=CaMa_dir+"/map/"+map+"/grdc_loc.txt"
with open(fname,"r") as f:
	lines=f.readlines()
outfile="../dat/grdc_"+map+".txt"
with open(outfile,"w") as fout:
    header="%7s;%34s;%42s;%7s;%7s;%7s;%7s\n"%("ID","River","Station","ix1","iy1","ix2","iy2")
    fout.write(header)
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print (line)
        num     = line[0].strip()
        basin   = line[1].strip()
        stream  = line[2].strip()
        ix1     = int(line[3])
        iy1     = int(line[4])
        ix2     = int(line[5])
        iy2     = int(line[6])
        staid   = int(num)
        org=grdc.grdc_dis(num,syear,eyear)
        org=np.array(org)
        # if np.sum(ma.masked_where(org!=-99.9,org))==0.0:
        if np.sum((org!=-9999.0)*1.0) < 365*thr_yr:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0))
            continue
        linew   = "%7s;%34s;%42s;%7d;%7d;%7d;%7d\n"%(num,basin,stream,ix1,iy1,ix2,iy2)
        print ("Obs >", 356*thr_yr,"days : ",num,basin,stream,np.sum((org!=-9999.0)*1.0))
        # print (linew)
        fout.write(linew)