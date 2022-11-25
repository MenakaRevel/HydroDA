#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from multiprocessing import Pool

# from gosh.params_real import ens_mem

# copy mean and std to HydroDA
#=============================
def copy(inputlist):
    input_name=inputlist[0]
    mapname=inputlist[1]
    tagout=inputlist[2]
    ens_char=input_name[-3::]
    runname="E2O" #pm.runname(pm.mode())
    # mean
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_sfcelv_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    oname="/cluster/data6/menaka/HydroDA/dat/mean_cal_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_sfcelv_cal_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    iname="/cluster/data6/menaka/ensemble_simulations/CaMa_out/"+input_name+"/sfcelv_mean"+tagout+".bin"

    print ("cp "+iname+" "+oname)
    os.system("cp "+iname+" "+oname)

    # std
    # oname="/cluster/data6/menaka/HydroDA/dat/std_sfcelv_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/std_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    oname="/cluster/data6/menaka/HydroDA/dat/std_cal_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/std_sfcelv_cal_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    iname="/cluster/data6/menaka/ensemble_simulations/CaMa_out/"+input_name+"/sfcelv_std"+tagout+".bin"

    print ("cp "+iname+" "+oname)
    os.system("cp "+iname+" "+oname)

    # print (ens_char, oname)
    return
#=============================
inputlist=[]
syear="%04d"%(2000) #2000)
eyear="%04d"%(2014) #2014) #2010)
tag=syear+"-"+eyear
# expname="AMZ049" #pm.expname()
expname="AMZCAL049"
runname="E2O" #pm.runname(pm.mode())
mapname="amz_06min" #pm.mapname()
ens_mem=49
for ens in np.arange(1,ens_mem+1): #pm.ens_mem(pm.mode())+1
    inputname="%s%s%03d"%(expname,runname,ens)
    inputlist.append([inputname,mapname,tag])
#==============
para=20
p=Pool(para)
p.map(copy,inputlist)
p.terminate()
#map(stats,inputlist)