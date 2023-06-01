#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
import os
import calendar
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
import math

# import CaMa-Flood variable reading using fortran
sys.path.append('../etc/')
from read_CMF import read_discharge, read_discharge_multi
#===============================================================================
# Experiment name
#===============================================================================
experiment="DIR_WSE_ISIMIP3a_SWOT_051" #"NOM_WSE_ERA5_CGLS_003"
#===============================================================================
# assim_out="../out/"+experiment
assim_out="/cluster/data7/menaka/HydroDA/out/"+experiment
print (assim_out)
#===============================================================================
# HydroDA related functions
sys.path.append(assim_out)
import params as pm
import read_grdc as grdc
import cal_stat as stat
#========================================
#====  functions for making figures  ====
#========================================
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    where ever the observed data contains nan
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    return data[:,0],data[:,1]
#========================================
def SWOT_day(yyyy,mm,dd):
  st_year,st_month,st_date=pm.starttime()
  start_time=datetime.date(st_year,st_month,st_date)
  this_time=datetime.date(int(yyyy),int(mm),int(dd))
  days=this_time-start_time
  days=days.days
  return days%21+1
#========================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#========================================
def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s) 
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)
#========================================
def KGE(s,o):
    """
	Kling Gupta Efficiency (Kling et al., 2012, http://dx.doi.org/10.1016/j.jhydrol.2012.01.011)
	input:
        s: simulated
        o: observed
    output:
        KGE: Kling Gupta Efficiency
    """
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    B = np.mean(s) / np.mean(o)
    y = (np.std(s) / np.mean(s)) / (np.std(o) / np.mean(o))
    r = np.corrcoef(o, s)[0,1]
    return 1 - np.sqrt((r - 1) ** 2 + (B - 1) ** 2 + (y - 1) ** 2)
#========================================
def correlation(s,o):
    """
    correlation coefficient
    input:
        s: simulated
        o: observed
    output:
        correlation: correlation coefficient
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    if s.size == 0:
        corr = 0.0 #np.NaN
    else:
        corr = np.corrcoef(o, s)[0,1]
        
    return corr
#========================================
def RMSE(s,o):
    """
    Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        RMSE: Root Mean Squre Error
    """
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    s,o = filter_nan(s,o)
    # return np.sqrt(np.mean((s-o)**2))
    return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))
#====================================================================
def read_dis(ix1, iy1, ix2, iy2, syear, eyear, indir):
    """
    Read CaMa-Flood discharge
    """
    dis = np.zeros( nbdays, 'f')
    # dis_max = np.zeros( nbyears, 'f')
    for year in range(syear, eyear):
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
        
        f = indir + '/outflw'+str(year)+'.bin'

        #outflw = np.fromfile(f, 'float32').reshape(-1,len(lat_global), len(lon_global))
        #if ix2 < 0 and iy2 < 0:
        #    tmp = outflw[:,iy1, ix1]
        #else:
        #    tmp = outflw[:,iy1, ix1] + outflw[:,iy2, ix2]

        #print tmp

        tmp = read_discharge( ix1+1, iy1+1, ix2+1, iy2+1, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        dis[s_days:e_days] = tmp

    return dis
#====================================================================
def read_dis_multi(ix1, iy1, ix2, iy2, syear, eyear, indir):
    #print ix1,iy1
    dis = np.zeros( (len(ix1), nbdays), 'f')
    dis_max = np.zeros( (len(ix1), nbyears), 'f')
    for year in range(syear, eyear+1):
        print (year)
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
        
        f = indir + '/outflw'+str(year)+'.bin'

        tmp = read_discharge_multi( ix1, iy1, ix2, iy2, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        dis[:,s_days:e_days] = tmp
        dis_max[:,year-syear] = np.nanmax(tmp, axis=1)

    return dis , dis_max
#====================================================================
mk_dir(assim_out+"/figures")
mk_dir(assim_out+"/figures/disgraph")
#----
fname=pm.CaMa_dir()+"/map/"+pm.mapname()+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
#----
syear,smonth,sdate=2001,1,1 #pm.starttime()
eyear,emonth,edate=2002,1,1 #pm.endtime()
#month=1
#date=1
start_dt=datetime.date(syear,smonth,sdate)
end_dt=datetime.date(eyear,emonth,edate)
size=60

start=0
last=(end_dt-start_dt).days
nbdays=int(last)
#last=365#int(argvs[1])
#if calendar.isleap(year):
#    last=366
#else:
#    last=365

ncpus=20
#last=89
N=int(last)
print ("days: ",N)
green2="greenyellow"
green ="green"
#colors = pc.color_assimd()
staid=[]
pname=[]
xlist=[]
ylist=[]
river=[]
#--
# rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZON","MEKONG","IRRAWADDY","VOLGA", "NIGER","YUKON","DANUBE"] #,"INDUS"] #["AMAZONAS"]#["CONGO"]#
# rivernames  = ["AMAZON"]
# rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZON","IRRAWADDY","VOLGA","NIGER","YUKON","DANUBE"] #,"INDUS"] #["AMAZONAS"]#["CONGO"]#
rivernames  = ["AMAZON", "MISSISSIPPI","MEKONG", "VOLGA", "COLORADO","MISSOURI","NIGER"]
#rivernames  = ["AMAZON"]
# rivernames  = ["COLORADO"]
# rivernames  = ["CHURCHILL"]
# rivernames = ["SAINT LAWRENCE","OHIO","CONNECTICUT","MISSOURI","MISSISSIPPI","COLORADO","CHURCHILL"]
# rivernames = grdc.grdc_river_name_v396()
for rivername in rivernames:
  grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername,fname=pm.CaMa_dir() + "/map/"+pm.mapname()+"/grdc_loc.txt")
  print (rivername, grdc_id,station_loc)
  river.append([rivername]*len(station_loc))
  staid.append(grdc_id)
  pname.append(station_loc)
  xlist.append(x_list)
  ylist.append(y_list)
#--
river=([flatten for inner in river for flatten in inner])
staid=([flatten for inner in staid for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])

print (len(pname), len(xlist))

pnum=len(pname)

org=[]
opn=[]
asm=[]

swt={}
for point in np.arange(pnum):
    swt[point] = []
# multiprocessing array
# result       = np.ctypeslib.as_ctypes(np.zeros((size, size)))
# shared_array = sharedctypes.RawArray(result._type_, result)
opn=np.ctypeslib.as_ctypes(np.zeros([N,pm.ens_mem(),pnum],np.float32))
shared_array_opn  = sharedctypes.RawArray(opn._type_, opn)
asm=np.ctypeslib.as_ctypes(np.zeros([N,pm.ens_mem(),pnum],np.float32))
shared_array_asm  = sharedctypes.RawArray(asm._type_, asm)
# for parallel calcualtion
inputlist=[]
for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        inputlist.append([yyyy,mm,dd,numch])
        # print (yyyy,mm,dd,numch)

def read_data(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    numch= inputlist[3]
    print (yyyy,mm,dd,numch)
    #--
    tmp_opn  = np.ctypeslib.as_array(shared_array_opn)
    tmp_asm  = np.ctypeslib.as_array(shared_array_asm)

    # year, mon, day
    year=int(yyyy)
    mon=int(mm)
    day=int(dd)
    num=int(numch)-1
    #--
    target_dt=datetime.date(year,mon,day)
    dt=(target_dt-start_dt).days
    # corrpted
    fname=assim_out+"/assim_out/outflw/open/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    print ("---- reading file ->", fname)
    #fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=assim_out+"/assim_out/outflw/assim/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    print ("---- reading file ->", fname)
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #-------------
    for point in np.arange(pnum):
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point],fname=pm.CaMa_dir() + "/map/"+pm.mapname()+"/grdc_loc.txt")
        if ix2 == -9999 or iy2 == -9999:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]
        else:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]+opnfile[iy2,ix2]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]+asmfile[iy2,ix2]
#--------
p   = Pool(ncpus)
res = p.map(read_data, inputlist)
opn = np.ctypeslib.as_array(shared_array_opn)
asm = np.ctypeslib.as_array(shared_array_asm)
p.terminate()

#-----------------------------------
# res = map(read_data, inputlist)
# opn = np.ctypeslib.as_array(shared_array_opn)
# asm = np.ctypeslib.as_array(shared_array_asm)

# for day in np.arange(start,last):
#     target_dt=start_dt+datetime.timedelta(days=day)
#     yyyy='%04d' % (target_dt.year)
#     mm='%02d' % (target_dt.month)
#     dd='%02d' % (target_dt.day)
#     print yyyy,mm,dd

# #    fname="../data/mesh_day%02d.bin"%(SWOT_day(yyyy,mm,dd))
# #    mesh_in=np.fromfile(fname,np.float32).reshape([640,1440])
# #    mesh=(mesh_in>=10)*(mesh_in<=60)
# #    meshP=mesh-1000*(mesh<0.1)

# #    fname="../sat/observation_day%02d.bin"%(SWOT_day(yyyy,mm,dd))
# #    meshP=np.fromfile(fname,np.int32).reshape([ny,nx])
# #    #meshP=(meshP>=1)*1


# #    # make org
# #    fname=assim_out+"/rivout/true/rivout"+yyyy+mm+dd+".bin"
# #    orgfile=np.fromfile(fname,np.float32).reshape([ny,nx])
# #
# #    org_frag=[]
# #    for point in np.arange(pnum):
# #        #print point
# #        xpoint=xlist[point]
# #        ypoint=ylist[point]
# #        org_frag.append(orgfile[ypoint,xpoint])
# #        #---SWOT--
# #        #if meshP[ypoint-40,xpoint] >= 1:
# #        if meshP[ypoint,xpoint] >= 1:
# #          if point not in swt.keys():
# #            swt[point] = [day]
# #          else:
# #            swt[point].append(day)
# #
# #    org.append(org_frag)

#     # make asm and opn
#     opn_ens=[]
#     asm_ens=[]
#     for num in np.arange(1,pm.ens_mem()+1):
#         numch='%03d'%num

#         fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
#         #fname="../CaMa_out/"+yyyy+mm+dd+"C"+numch+"/rivout"+yyyy+".bin"
#         opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])

#         fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
#         #fname="../CaMa_out/"+yyyy+mm+dd+"A"+numch+"/rivout"+yyyy+".bin"
#         asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])

#         opn_frag=[]
#         asm_frag=[]
#         for point in np.arange(pnum):
#             ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point])
#             xpoint=xlist[point]
#             ypoint=ylist[point]
#             if ix2 == -9999 or iy2 == -9999:
#                 opn_frag.append(opnfile[iy1,ix1])
#                 asm_frag.append(asmfile[iy1,ix1])
#             else:
#                 opn_frag.append(opnfile[iy1,ix1]+opnfile[iy2,ix2])
#                 asm_frag.append(asmfile[iy1,ix1]+asmfile[iy2,ix2])
#             #opn_frag.append(opnfile[ypoint,xpoint])
#             #asm_frag.append(asmfile[ypoint,xpoint])

#         opn_ens.append(opn_frag)
#         asm_ens.append(asm_frag)


#     opn.append(opn_ens)
#     asm.append(asm_ens)

# #org=np.array(org,dtype=np.float32)
# opn=np.array(opn,dtype=np.float32)
# asm=np.array(asm,dtype=np.float32)

##--
#print np.shape(org),org.dtype
#org.tofile("org.bin")
#print np.shape(opn)
#opn.tofile("opn.bin")
#print np.shape(asm)
#asm.tofile("asm.bin")
##--
#def save_txt(data,name):
#  data=data.flatten()
#  f=open(name,"w")
#  for i in data:
#    line="%10.4f\n"%(i)
#    f.write(line)
#  f.close()
#  return 0
#
#save_txt(org,"org.txt")
#save_txt(opn,"opn.txt")
#save_txt(asm,"asm.txt")

#for point in np.arange(pnum):
def make_fig(point):
    plt.close()
    #labels=["GRDC","corrupted","assimilated"]
    obstype=pm.obs_name()
    if obstype=="SWOT":
        exptype="virtual"
        labels=["true","simulated","assimilated"]
    else:
        exptype="real"
        labels=["GRDC","simulated","assimilated"]
    #
    #print org[:,point]
    #for i in np.arange(start,last):
        #print opn[i,:,point]
        #print asm[i,:,point]

#    plt.plot(np.arange(start,last),org[:,point],label="true",color=colors["true"],linewidth=0.7)
#
#    for num in np.arange(0,pm.ens_mem()):
#        plt.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color=colors["corrupted"],linewidth=0.3,alpha=0.5)
#        plt.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color=colors["assimilated"],linewidth=0.3,alpha=0.5)
#
#    plt.ylim(ymin=0)
    if exptype=="virtual":
        # org=read_dis()
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point])
        # indir = "/work/a04/julien/CaMa-Flood_v4/out/coupled-model2"
        # indir = "/cluster/data6/menaka/CaMa-H08/out/obs_org"
        # indir = "/cluster/data6/menaka/CaMa-H08/out/obs_rivhgt"
        # indir = "/cluster/data6/menaka/CaMa-H08/out/obs_rivwth"
        indir = "/cluster/data6/menaka/CaMa-H08/out/obs_rivman"
        # indir = "/cluster/data6/menaka/CaMa-H08/out/obs_fldhgt"
        # indir = "/cluster/data6/menaka/CaMa-H08/out/obs_corr_all"
        # indir = pm.obs_dir()
        org=read_dis(ix1, iy1, ix2, iy2, syear, eyear, indir)
        # org=grdc.grdc_dis(staid[point],syear,eyear-1)
        org=np.array(org)
    else:
        org=grdc.grdc_dis(staid[point],syear,eyear-1) #,smon=1,emon=1,sday=1,eday=31)
        org=np.array(org)
    #------------------------
    # Making Figure
    fig, ax1 = plt.subplots()
    lines=[ax1.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="#34495e",linewidth=3.0,zorder=101)[0]] #,marker = "o",markevery=swt[point])
#    ax1.plot(np.arange(start,last),hgt[:,point],label="true",color="gray",linewidth=0.7,linestyle="--",zorder=101)
#    plt.plot(np.arange(start,last),org[:,point],label="true",color="black",linewidth=0.7)
    for num in np.arange(0,pm.ens_mem()):
        ax1.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color="blue",linewidth=0.1,alpha=0.1,zorder=102)
        ax1.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color="red",linewidth=0.1,alpha=0.1,zorder=103)
#        plt.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color="blue",linewidth=0.3,alpha=0.5)
#        plt.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color="red",linewidth=0.3,alpha=0.5)
    # draw mean of ensembles
    lines.append(ax1.plot(np.arange(start,last),np.mean(ma.masked_less(opn[:,:,point],0.0),axis=1),label="corrupted",color="#4dc7ec",linewidth=1.0,alpha=1,zorder=104)[0])
    lines.append(ax1.plot(np.arange(start,last),np.mean(ma.masked_less(asm[:,:,point],0.0),axis=1),label="assimilated",color="#ff8021",linewidth=1.0,alpha=1,zorder=106)[0])
    #    plt.ylim(ymin=)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('discharge (m$^3$/s)', color='k')
    ax1.set_xlim(xmin=0,xmax=last+1)
    ax1.tick_params('y', colors='k')
    # scentific notaion
    ax1.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax1.yaxis.major.formatter._useMathText=True 
    #
    #xxlist=np.linspace(0,N,(eyear-syear)+1)
    #xlab=np.arange(syear,eyear+1,1)
    #xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    if eyear-syear > 5:
        dtt=2
        dt=int(math.ceil(((eyear-syear)+1)/2.0))
    elif eyear-syear > 10:
        dtt=5
        dt=int(math.ceil(((eyear-syear)+1)/5.0))
    else:
        dtt=1
        dt=(eyear-syear)+1
    xxlist=np.linspace(0,N,dt,endpoint=True)
    #xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    xxlab=np.arange(syear,eyear+1,dtt)
    ax1.set_xticks(xxlist)
    ax1.set_xticklabels(xxlab,fontsize=10)
    # Nash-Sutcllf calcuation
    NS1=NS(np.mean(asm[:,:,point],axis=1),org)
    NS2=NS(np.mean(opn[:,:,point],axis=1),org)
    KGE1=KGE(np.mean(asm[:,:,point],axis=1),org)
    KGE2=KGE(np.mean(opn[:,:,point],axis=1),org)
    COR1=correlation(np.mean(asm[:,:,point],axis=1),org)
    COR2=correlation(np.mean(opn[:,:,point],axis=1),org)
    RSE1=RMSE(np.mean(asm[:,:,point],axis=1),org)
    RSE2=RMSE(np.mean(opn[:,:,point],axis=1),org)
    #1-((np.sum((org[:ed,point]-org_Q)**2))/(np.sum((org_Q-np.mean(org_Q))**2)))
    #print point,NS1,NS2
    Nash1="NS (assim):%4.2f"%(NS1)
    Nash2="NS (open):%4.2f"%(NS2)
    kgeh1="KGE(assim):%4.2f"%(KGE1)
    kgeh2="KGE(open):%4.2f"%(KGE2)
    corr1="$r$(assim):%4.2f"%(COR1)
    corr2="$r$(open):%4.2f"%(COR2)
    rmse1="RMSE: %4.2f"%(RSE1)
    rmse2="RMSE: %4.2f"%(RSE2)
    # #
    # ax1.text(0.02,0.95,Nash1,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    # ax1.text(0.02,0.85,Nash2,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    # ax1.text(0.42,0.95,kgeh1,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    # ax1.text(0.42,0.85,kgeh2,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    # ax1.text(0.02,0.75,corr1,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    # ax1.text(0.02,0.55,corr2,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    #
    ax1.text(0.02,0.95,Nash1,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.02,0.90,Nash2,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.02,0.80,kgeh1,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.02,0.75,kgeh2,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.42,0.95,corr1,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.42,0.90,corr2,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.42,0.80,rmse1,ha="left",va="center",transform=ax1.transAxes,fontsize=8)
    ax1.text(0.42,0.75,rmse2,ha="left",va="center",transform=ax1.transAxes,fontsize=8)

#    # twin axis
#    ax2 = ax1.twinx()
#    #aiv = stat.AI(asm[:,:,point],opn[:,:,point],org[:,point])
#    #aivn= stat.AI_new(asm[:,:,point],opn[:,:,point],org[:,point])
#    aivn,error = stat.AI_new(asm[:,:,point],opn[:,:,point],org[:,point])
#    #print aivn
#    #--
#    pBias,pB_c = stat.pBias(asm[:,:,point],opn[:,:,point],org[:,point])
#    ai_mean =np.mean(ma.masked_less_equal(aivn,0.0))# np.nanmean(ma.masked_less_equal(aivn,0.0))
#    # RMSE
#    RootMSE=stat.RMSE(asm[:,:,point],org[:,point])
#    # rRMSE
#    rRootMSE=stat.rRMSE(asm[:,:,point],org[:,point])
#    # NRMSE
#    NRootMSE=stat.NRMSE(asm[:,:,point],org[:,point])
#    # VE
#    VolEff=stat.VE(asm[:,:,point],org[:,point])
#    # NSE
#    NSEc,NSa,NSc=stat.NSE(asm[:,:,point],opn[:,:,point],org[:,point])
#    # PDRI PTRI
#    PDRI,PTRI=stat.PRI(asm[:,:,point],opn[:,:,point],org[:,point])
#    #---
#    EnsSprd=stat.EnsSpr(asm[:,:,point],opn[:,:,point],org[:,point])
#    EnsSpr_mean=np.mean(ma.masked_less_equal(EnsSprd,0.0))
#    #---
#    AssimQlt=stat.AQ(PDRI,PTRI,ai_mean,EnsSpr_mean)
#    #---
#    mai = "meanAI:%1.2f"%(ai_mean)
#    pB  = "pBIAS:%1.1f%%"%(pBias)
#    rms = "RMSE:%8.2f"%(RootMSE)
#    rrms= "rRMSE:%3.2f"%(rRootMSE)
#    nrms= "NRMSE:%3.2f"%(NRootMSE)
#    veff= "VE:%3.2f"%(VolEff)
#    nsec= "NSE:%3.2f"%(NSEc)
#    pdr = "PDRI:%3.2f"%(PDRI)
#    ptr = "PTRI:%3.2f"%(PTRI)
#    ESrd= "EnsSpr:%3.2f"%(EnsSpr_mean)
#    AssQ= "AQ:%3.2f"%(AssimQlt)
#    NSac= "NSEa:%3.2f, NSEc:%3.2f"%(NSa,NSc)
#    print mai , pB#, "pB_c", pB_c, "%"
#    #pB  = pB + r'{1.1f}\%'.format(pBias)
#    ax1.text(0.02,0.95,mai,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.85,nsec,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.75,pB,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.65,rms,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.55,rrms,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.45,nrms,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.35,veff,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.25,pdr,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.15,ptr,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    #ax1.text(0.02,0.05,AssQ,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.05,NSac,ha="left",va="center",transform=ax1.transAxes,fontsize=10)

    #---
#    pBA_p=np.sum(ma.masked_greater((np.mean(asm[:,:,point],axis=1)-org[:,point]),0.0).filled(0.0))
#    pBA_n=np.sum(ma.masked_less((np.mean(asm[:,:,point],axis=1)-org[:,point]),0.0).filled(0.0))
#    #---
#    pBC_p=np.sum(ma.masked_greater((np.mean(opn[:,:,point],axis=1)-org[:,point]),0.0).filled(0.0))
#    pBC_n=np.sum(ma.masked_less((np.mean(opn[:,:,point],axis=1)-org[:,point]),0.0).filled(0.0))
#    ax1.text(0.80,0.95,"pBa_p:%1.1f%%"%(pBA_p),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.80,0.85,"pBa_n:%1.1f%%"%(pBA_n),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.80,0.75,"pBc_p:%1.1f%%"%(pBC_p),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.80,0.65,"pBc_n:%1.1f%%"%(pBC_n),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.80,0.55,"pBa_p/pBa_n:%1.1f"%(abs(pBA_p/pBA_n+1.0e-20)),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.80,0.45,"pBc_p/pBc_n:%1.1f"%(abs(pBC_p/pBC_n+1.0e-20)),ha="left",va="center",transform=ax1.transAxes,fontsize=10)


#    ax2.plot(np.arange(start,last),aiv,color="green",zorder=104,marker = "o",alpha=0.3, markevery=swt[point])
#    ax2.plot(np.arange(start,last),ma.masked_less_equal(aivn,1.0e-20),color="green",zorder=104,marker = "o", alpha=0.5,linewidth=0.5,markevery=swt[point])
#    ax2.plot(np.arange(start,last),ma.masked_where(error<=0.1,aivn),color=green,marker = "o",markeredgecolor =green, alpha=1,linewidth=0.5,markevery=swt[point],markersize=3,zorder=100)#,label="AI")
#    ax2.plot(np.arange(start,last),ma.masked_where(error>0.1,aivn),color=green2,marker = "o",markeredgecolor =green2, alpha=1,linewidth=0.5,markevery=swt[point],markersize=3,zorder=100)#,label="AI")
#    ax2.set_ylabel('AI', color='green')
#    ax2.tick_params('y', colors='green')
#    ax2.set_ylim(ymin=0.,ymax=1.)
#    ax2.set_xlim(xmin=0,xmax=last+1)
#    print swt[point]
    plt.legend(lines,labels,ncol=1,loc='upper right') #, bbox_to_anchor=(1.0, 1.0),transform=ax1.transAxes)
    station_loc_list=pname[point].split("/")
    print (station_loc_list)
    station_name="".join(station_loc_list[0].split())
    # station_name="-".join(station_loc_list) 
    print ('--- saving figure',river[point]+"-"+station_name+".png")
    plt.savefig(assim_out+"/figures/disgraph/"+river[point]+"-"+station_name+".png",dpi=500)
    return 0



#p=Pool(12)
#p.map(make_fig,np.arange(pnum))
#p.terminate()


#plt.clf()
#fig = plt.figure()
#figcbr = plt.figure(figsize=(8,0.6))
#axcbr  = fig.add_axes([0,0.5,1.0,0.6])
#plt.plot([],[],label="true",color=colors["true"],linewidth=1)
#plt.plot([],[],label="corrupted",color=colors["corrupted"],linewidth=1,alpha=0.5)
#plt.plot([],[],label="assimilated",color=colors["assimilated"],linewidth=1,alpha=0.5)
#plt.legend(loc=10,ncol=3)
#plt.savefig("../assim_out/fig/disgraph/legend.png",dpi=500)




para_flag=1
# para_flag=0
#--
if para_flag==1:
    p=Pool(ncpus)
    p.map(make_fig,np.arange(pnum))
    p.terminate()
else:
    map(make_fig,np.arange(pnum))
