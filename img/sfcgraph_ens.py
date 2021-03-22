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

os.system("ln -sf ../gosh/params.py params.py")
#sys.path.append('../assim_out/')
import params as pm
import read_grdc as grdc
import read_hydroweb as hweb
import cal_stat as stat
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

#argvs = sys.argv

#experiment="E2O_HydroWeb22"
# experiment="VIC_BC_HydroWeb11"
# experiment="test_wse"
# experiment="DIR_WSE_E2O_HWEB_001"
experiment="ANO_WSE_E2O_HWEB_001"
# experiment="NOM_WSE_E2O_HWEB_001"
conflag=2
#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
assim_out=pm.DA_dir()+"/out/"+experiment
print (assim_out)

#assim_out="assim_out_E2O_womc"
#assim_out="assim_out"
#assim_out="assim_out_E2O_wmc"
#assim_out="assim_out_E2O_womc_0"
#assim_out="assim_out_ECMWF_womc_baised_0"
#assim_out="assim_out_ECMWF_womc_baised"
#assim_out="assim_out_ECMWF_womc_baised_if"
#assim_out="assim_out_ECMWF_womc_baised_0.75"
#assim_out="assim_out_ECMWF_womc_baised_0.80"
#assim_out="assim_out"
#assim_out="assim_out_biased_womc"
#assim_out="assim_out_biased_wmc"


#os.system("mkdir ../assim_out/img")
#os.system("mkdir ../assim_out/img/sfcelv")
#----
def SWOT_day(yyyy,mm,dd):
  st_year,st_month,st_date=pm.starttime()
  start_time=datetime.date(st_year,st_month,st_date)
  this_time=datetime.date(int(yyyy),int(mm),int(dd))
  days=this_time-start_time
  days=days.days
  return days%21+1
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#----
mk_dir(assim_out+"/figures")
mk_dir(assim_out+"/figures/sfcelv")
#----
fname=pm.CaMa_dir()+"/map/"+pm.mapname()+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
#---
###year,month,date=pm.starttime()
####month=1
####date=1
###start_dt=datetime.date(year,month,date)
###size=60
###
###start=0
####last=int(argvs[1])
###last=365#int(argvs[1])
###if calendar.isleap(year):
###    last=366
###else:
###    last=365
syear,smonth,sdate=2002,1,1 #pm.starttime()#2004#1991
eyear,emonth,edate=pm.endtime() #2005,1,1 #
#month=1
#date=1
start_dt=datetime.date(syear,smonth,sdate)
end_dt=datetime.date(eyear,emonth,edate)
size=60

start=0
last=(end_dt-start_dt).days

N=int(last)
green2="greenyellow"
#--------------
nextxy = pm.CaMa_dir()+"/map/"+pm.mapname()+"/nextxy.bin"
rivwth = pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivwth_gwdlr.bin"
rivhgt = pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivhgt.bin"
rivlen = pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivlen.bin"
elevtn = pm.CaMa_dir()+"/map/"+pm.mapname()+"/elevtn.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
#----
# mean
#mean_sfcelv = pm.DA_dir()+"/dat/mean_sfcelv_1958-2013.bin"
#mean_sfcelv = pm.DA_dir()+"/dat/mean_sfcelv_E2O_amz_06min_1980-2014.bin"
#mean_sfcelv = pm.DA_dir()+"/dat/mean_sfcelv_VIC_BC_1980-2014.bin"
#mean_sfcelv = pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1980-2014.bin"
#mean_sfcelv = pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin"
mean_sfcelv = assim_out+"/assim_out/mean_sfcelv/mean_sfcelv.bin"
mean_sfcelv = np.fromfile(mean_sfcelv,np.float32).reshape(ny,nx)
# std
#std_sfcelv = pm.DA_dir()+"/dat/std_sfcelv_1958-2013.bin"
#std_sfcelv = pm.DA_dir()+"/dat/std_sfcelv_E2O_amz_06min_1980-2014.bin"
#std_sfcelv = pm.DA_dir()+"/dat/std_sfcelv_VIC_BC_1980-2014.bin"
#std_sfcelv = pm.DA_dir()+"/dat/std_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1980-2014.bin"
#std_sfcelv = pm.DA_dir()+"/dat/std_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin"
std_sfcelv = assim_out+"/assim_out/mean_sfcelv/std_sfcelv.bin"
std_sfcelv = np.fromfile(std_sfcelv,np.float32).reshape(ny,nx)
#- mean obs HydroWeb
mean_obs = pm.HydroWeb_dir()+"/bin/HydroWeb_mean.bin"
mean_obs = np.fromfile(mean_obs,np.float32).reshape(ny,nx)
# std obs HydroWeb
std_obs = pm.HydroWeb_dir()+"/bin/HydroWeb_std.bin"
std_obs = np.fromfile(std_obs,np.float32).reshape(ny,nx)
#-------
# # mean & std from previous year
# mean_obss=np.zeros([pm.ens_mem(),ny,nx])
# std_obss=np.zeros([pm.ens_mem(),ny,nx])
# for num in np.arange(1,int(pm.ens_mem())+1):
#     numch='%03d'%num
#     fname=assim_out+"/assim_out/mean_sfcelv/meansfcelvC"+numch+".bin"
#     mean_corr=np.fromfile(fname,np.float32).reshape([ny,nx])
#     mean_obss[num-1]=mean_corr
#     fname=assim_out+"/assim_out/mean_sfcelv/stdsfcelvC"+numch+".bin"
#     std_corr=np.fromfile(fname,np.float32).reshape([ny,nx])
#     std_obss[num-1]=std_corr
# mean_sfcelv=np.mean(mean_obss,axis=0)
# std_sfcelv=np.std(std_obss,axis=0)
###------
pname=[]
xlist=[]
ylist=[]
river=[]
EGM08=[]
EGM96=[]
#--
#rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS","MEKONG","IRRAWADDY","VOLGA", "NIGER","YUKON","DANUBE"] #,"INDUS"] #["AMAZONAS"]#["CONGO"]#
rivernames  = ["AMAZONAS"]
for rivername in rivernames:
  path = assim_out+"/figures/sfcelv/%s"%(rivername)
  #print path
  #mk_dir(path)
  #station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
  station_loc,x_list,y_list,egm08,egm96 =hweb.get_hydroweb_loc(rivername,pm.mapname())
  #print rivername, station_loc
  river.append([rivername]*len(station_loc))
  pname.append(station_loc)
  xlist.append(x_list)
  ylist.append(y_list)
  EGM08.append(egm08)
  EGM96.append(egm96)

##rivernames = grdc.grdc_river_name()
##for rivername in rivernames: 
##  station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
##  #print rivername, station_loc
##  for station in station_loc:
##    gid=grdc.get_id(station)
##    if gid== -9999:
##      continue 
##    path = "../"+assim_out+"/img/sfcelv/%s"%(rivername)
##    #print path
##    mk_dir(path)
##    ix, iy = grdc.get_loc_v394(gid)
##    print station,gid, ix ,iy
##    river.append([rivername])
##    pname.append([station])
##    xlist.append([ix])
##    ylist.append([iy])
###--
##  if rivername=="LENA":
##      river.append([rivername]*3)
##      pname.append(["L1","L2","L3"])
##      xlist.append([1233,1218,1206])
##      ylist.append([  71,  98, 117])
##  if rivername=="NIGER":
##      river.append([rivername]*6)
##      pname.append(["N1","N2","N3","N4","N5","N6"])
##      xlist.append([744,744,732,712,704,700])
##      ylist.append([342,324,310,292,295,303])
##  if rivername=="AMAZONAS":
##      river.append([rivername]*4)
##      pname.append(["B","E","F","G"])
##      xlist.append([515,447,464,420])
##      ylist.append([364,366,416,367])
##  if rivername=="MEKONG":
##      river.append([rivername]*3)
##      pname.append(["Me_D","Me_M","Me_U"])
##      xlist.append([1143,1139,1127])
##      ylist.append([ 319, 293, 282])
##  if rivername=="MISSISSIPPI":
##      river.append([rivername]*4)
##      pname.append(["Mi_D","Mi_M","Mi_U","Mi_U2"]) #最後はミズーリ川
##      xlist.append([361,362,345,308])
##      ylist.append([241,214,137,167])
##  if rivername=="OB":
##      river.append([rivername]*3)
##      pname.append(["O_D","O_M","O_U"])
##      xlist.append([995,996,1048])
##      ylist.append([ 92,121, 159])
###  if rivername=="CONGO":
###      pname.append(["C_D","C_M","C_U"])#,"C1","C2","C3","C4","C5","C6","C7","C8"])
###      xlist.append([772,813,834])#,769,784,789,809,821,824,802,794])
###      ylist.append([383,353,397])#,383,372,363,343,359,376,379,342])
###      river.append#([rivername,rivername,rivername])#,rivername,rivername,rivername,rivername,rivername,rivername,rivername,rivername])
##  if rivername=="INDUS":
##      river.append([rivername]*3)
##      pname.append(["I_D","I_Sub","I_M"])
##      xlist.append([992,995,1003])
##      ylist.append([251,252,233])
###
###  if rivername=="CONGO":
###    river.append([rivername]*10)
###    pname.append(["C4","C1","C10","C9","C7","C5","C3","C2","C8","C6"])
###    xlist.append([813,834,772,784,789,809,821,824,802,794])
###    ylist.append([353,397,383,372,363,343,359,376,379,342])
###
##  if rivername=="CONGO":
##      river.append([rivername]*6)
##      pname.append(["C2","C1","C6","C5","C3","C4"])
##      xlist.append([813,834,772,784,794,799])
##      ylist.append([353,397,383,372,342,375])

river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
EGM08=([flatten for inner in EGM08 for flatten in inner])
EGM96=([flatten for inner in EGM96 for flatten in inner])
pnum=len(pname)
#print len(river),pnum,pname,river
org=[]
opn=[]
asm=[]
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
        #print yyyy,mm,dd,numch

def read_data(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    numch= inputlist[3]
    #print (yyyy,mm,dd,numch)
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
    fname=assim_out+"/assim_out/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
    #fname=assim_out+"/assim_out/outflw/open/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=assim_out+"/assim_out/ens_xa/assim/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
    #fname=assim_out+"/assim_out/outflw/assim/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #print (yyyy,mm,dd,numch,dt)
    #-------------
    for point in np.arange(pnum):
        ix1=xlist[point]
        iy1=ylist[point]
        #print (ix1,iy1)
        tmp_opn[dt,num,point]=opnfile[iy1,ix1]
        tmp_asm[dt,num,point]=asmfile[iy1,ix1]
#--------
#print (pname[0])
p   = Pool(20)
res = p.map(read_data, inputlist)
#map(read_data, inputlist)
opn = np.ctypeslib.as_array(shared_array_opn)
asm = np.ctypeslib.as_array(shared_array_asm)
p.terminate()
#############
# hgt=[]
# bathy=[]
# ele=[]
# m_sf=[]
# em_sf=[]
# swt={}
# for point in np.arange(pnum):
#     swt[point] = []

# for day in np.arange(start,last):
#     target_dt=start_dt+datetime.timedelta(days=day)
#     yyyy='%04d' % (target_dt.year)
#     mm='%02d' % (target_dt.month)
#     dd='%02d' % (target_dt.day)
#     print yyyy,mm,dd

# #    fname="../sat/mesh_day%02d.bin"%(SWOT_day(yyyy,mm,dd))
# #    mesh_in=np.fromfile(fname,np.float32).reshape([640,1440])
# #    mesh=(mesh_in>=10)*(mesh_in<=60)
# #    meshP=mesh-1000*(mesh<0.1)
# #
# #    # make org
# #    fname=assim_out+"/xa_m/true/"+yyyy+mm+dd+"_xam.bin"
# #    #fname="../CaMa_out/"+yyyy+mm+dd+"T000/sfcelv"+yyyy+".bin"
# #    orgfile=np.fromfile(fname,np.float32).reshape([720,1440])
# #
# #    fname=pm.CaMa_dir()+"/map/glb_15min/rivhgt.bin"
# #    bathyfile=np.fromfile(fname,np.float32).reshape([720,1440])
# #
# #    fname=assim_out+"/mean_sfcelv/meansfcelvT000.bin"
# #    mean_true=np.fromfile(fname,np.float32).reshape([720,1440])
# #
# #    org_frag=[]
# #    bathy_frag=[]
# #    ele_frag=[]
# #    m_sf_frag=[]
# #    for point in np.arange(pnum):
# #        xpoint=xlist[point]
# #        ypoint=ylist[point]
# #        org_frag.append(orgfile[ypoint,xpoint])
# #        bathy_frag.append(elevtn[ypoint,xpoint] -bathyfile[ypoint,xpoint])
# #        ele_frag.append(elevtn[ypoint,xpoint])
# #        m_sf_frag.append(mean_true[ypoint,xpoint])
# #
# #        #---SWOT--
# #        if meshP[ypoint-40,xpoint] >= 1:
# #          if point not in swt.keys():
# #            swt[point] = [day]
# #          else:
# #            swt[point].append(day)
# #
# #    org.append(org_frag)
# #    bathy.append(bathy_frag)
# #    ele.append(ele_frag)
# #    m_sf.append(m_sf_frag)

#     # make asm and opn
#     opn_ens=[]
#     asm_ens=[]
#     hgt_ens=[]
#     em_sf_ens=[]
#     for num in np.arange(1,int(pm.ens_mem())+1):
#         numch='%03d'%num

#         fname=assim_out+"/assim_out/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
#         #fname="../CaMa_out/"+yyyy+mm+dd+"C"+numch+"/sfcelv"+yyyy+".bin"
#         #fname="../"+assim_out+"/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
#         opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])

#         fname=assim_out+"/assim_out/ens_xa/assim/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
#         asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])

# #        fname="../assim_out/rivhgt/assim/rivhgt"+yyyy+mm+dd+"_"+numch+"A.bin"
#         #fname="../CaMa_out/"+yyyy+mm+dd+"A"+numch+"/sfcelv"+yyyy+".bin"
#         rhgtfile=np.fromfile(fname,np.float32).reshape([ny,nx])

#         #fname=assim_out+"/assim_out/mean_sfcelv/meansfcelvC"+numch+".bin"
#         mean_corr=np.fromfile(fname,np.float32).reshape([ny,nx])

#         opn_frag=[]
#         asm_frag=[]
#         hgt_frag=[]
#         em_sf_frag=[]
#         for point in np.arange(pnum):
#             xpoint=xlist[point]
#             ypoint=ylist[point]
#             opn_frag.append(opnfile[ypoint,xpoint])
#             asm_frag.append(asmfile[ypoint,xpoint])
#             hgt_frag.append(rhgtfile[ypoint,xpoint])
#             #print asmfile[ypoint,xpoint],elevtn[ypoint,xpoint] - rhgtfile[ypoint,xpoint],rhgtfile[ypoint,xpoint]
#             em_sf_frag.append(mean_corr[ypoint,xpoint])
#         opn_ens.append(opn_frag)
#         asm_ens.append(asm_frag)
#         hgt_ens.append(hgt_frag)
#         em_sf_ens.append(em_sf_frag)

#     opn.append(opn_ens)
#     asm.append(asm_ens)
#     hgt.append(hgt_ens)
#     em_sf.append(em_sf_ens)
# #-----
# #org=np.array(org)
# #bathy=np.array(bathy)
# #ele=np.array(ele)
# opn=np.array(opn)
# asm=np.array(asm)
# #hgt=np.array(hgt)
# #m_sf=np.array(m_sf)
# em_sf=np.array(em_sf)

# #print np.shape(org),org.dtype
# #org.tofile("org.bin")
# #print np.shape(opn)
# #opn.tofile("opn.bin")
# #print np.shape(asm)
# #asm.tofile("asm.bin")
#======================
#--
# def save_txt(data,name):
#   data=data.flatten()
#   f=open(name,"w")
#   for i in data:
#     line="%10.4f\n"%(i)
#     f.write(line)
#   f.close()
#   return 0

#save_txt(org,"org.txt")
#save_txt(opn,"opn.txt")
#save_txt(asm,"asm.txt")

print 'making figure'
#for point in np.arange(pnum):
def make_fig(point):
    plt.close()
    labels=["HydroWeb","corrupted","assimilated"]
#    print org[:,point]
#    print "----------------"
#    print np.mean(asm[:,:,point],axis=1)
#    for i in np.arange(start,last):
#        print opn[i,:,point]
#        print asm[i,:,point]

    fig, ax1 = plt.subplots()
    #ax1.plot(np.arange(start,last),org[:,point],label="true",color="black",linewidth=0.7,zorder=101,marker = "o",markevery=swt[point])
    time,org=hweb.HydroWeb_WSE(pname[point],syear,eyear)
    data=np.array(org)+np.array(EGM08[point])-np.array(EGM96[point])
    #alti=hweb.altimetry(pname[point]) - EGM08[point] + EGM96[point]
    #data=data-alti+elevtn[ylist[point],xlist[point]]
    #data=((data-np.mean(data))/np.std(data))*std_sfcelv[ylist[point],xlist[point]]+mean_sfcelv[ylist[point],xlist[point]]
    #data=(data-np.mean(data))+mean_sfcelv[ylist[point],xlist[point]]
    #---------
    #data=data-mean_obs[ylist[point],xlist[point]]+mean_sfcelv[ylist[point],xlist[point]]
    if conflag==1:
        data0=data
    elif conflag==2:
        data0=(data-mean_obs[ylist[point],xlist[point]])+mean_sfcelv[ylist[point],xlist[point]]
    elif conflag==3:
        data0=((data-mean_obs[ylist[point],xlist[point]])/(std_obs[ylist[point],xlist[point]]+1.0e-20))*std_sfcelv[ylist[point],xlist[point]]+mean_sfcelv[ylist[point],xlist[point]]
    lines=[ax1.plot(time,data0,label="obs",marker="o",color="#34495e",linewidth=0.0,zorder=101)[0]]
#    ax1.plot(np.arange(start,last),org[:,point],label="true",color="black",linewidth=0.7,zorder=101)
#    ax1.plot(np.arange(start,last),m_sf[:,point],label="mean sfcelv",color="black",linewidth=0.7,linestyle="--",zorder=107)
#    plt.plot(np.arange(start,last),org[:,point],label="true",color="black",linewidth=0.7)

    for num in np.arange(0,int(pm.ens_mem())):
        ax1.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color="#4dc7ec",linewidth=0.2,alpha=0.5,zorder=102)
        ax1.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color="#ff8021",linewidth=0.1,alpha=0.5,zorder=103)
#        ax1.plot(np.arange(start,last),em_sf[:,num,point],label="mean sfcelv",color="blue",linewidth=0.3,linestyle="--",alpha=0.5,zorder=103)
#        plt.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color="blue",linewidth=0.3,alpha=0.5)
#        plt.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color="red",linewidth=0.3,alpha=0.5)
    lines.append(ax1.plot(np.arange(start,last),np.mean(opn[:,:,point],axis=1),label="corrupted",color="#4dc7ec",linewidth=0.8,alpha=0.8,zorder=102)[0])
    lines.append(ax1.plot(np.arange(start,last),np.mean(asm[:,:,point],axis=1),label="assimilated",color="#ff8021",linewidth=0.8,alpha=0.8,zorder=103)[0])
#    ax1.plot(np.arange(start,last),np.mean(em_sf[:,:,point],axis=1),label="mean sfelv",color="blue",linewidth=0.5,linestyle="--",alpha=0.5,zorder=103)
#    plt.ylim(ymin=)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('WSE (m)', color='k')
    #ax1.set_ylim(ymin=0,ymax=250.)
    ax1.set_xlim(xmin=0,xmax=last+1)
    ax1.tick_params('y', colors='k')
    # xxlist=np.linspace(15,N-15,int(N/30))
    # xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    #ax1.set_xticks(xxlist)
    #ax1.set_xticklabels(xxlab,fontsize=10)
    #--for bug fixing
    corr_mean=np.mean(np.mean(opn[:,:,point],axis=1),axis=0)
    obs_mean=mean_obs[ylist[point],xlist[point]]
    sfc_mean=mean_sfcelv[ylist[point],xlist[point]]
    outtext="courrpted mean: %6.2f"%(corr_mean)
    ax1.text(0.02,0.9,outtext,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    outtext="used(norm) mean: %6.2f"%(sfc_mean)
    ax1.text(0.02,0.8,outtext,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    outtext="observation mean: %6.2f"%(obs_mean)
    ax1.text(0.02,0.7,outtext,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    # xlable in years
    if eyear-syear > 5:
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

#    ax2 = ax1.twinx()
#    aiv = stat.AI(asm[:,:,point],opn[:,:,point],org[:,point])
#    aivn,error= stat.AI_new(asm[:,:,point],opn[:,:,point],org[:,point])#[0]
#    fmax = np.finfo(np.float32).max
#    aivn[aivn>fmax]=0
#    aivn[aivn<-fmax]=0
#    print np.mean(aivn)
#    #--
#    pBias,pB_c = stat.pBias(asm[:,:,point],opn[:,:,point],org[:,point])
#    ai_mean = np.mean(ma.masked_less_equal(aivn,0.0)) #np.nanmean(ma.masked_less_equal(aivn,0.0))
#    mai = "meanAI:%1.2f"%(ai_mean)
#    pB  = "pBIAS:%1.1f%%"%(pBias)
#    print  mai # pB, "pB_c", pB_c, "%"
#    #pB  = pB + r'{1.1f}\%'.format(pBias)
#    ax1.text(0.02,0.9,mai,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#    ax1.text(0.02,0.8,pB,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
##    ax2.plot(np.arange(start,last),aiv,color="green",zorder=104,marker = "o",alpha=0.3, markevery=swt[point])
#    ax2.plot(np.arange(start,last),ma.masked_less_equal(aivn,1.0e-20),color="green",zorder=104,marker = "o", alpha=0.5,linewidth=0.5,markevery=swt[point])
#    ax2.set_ylabel('AI', color='green')
#    ax2.tick_params('y', colors='green')
#    ax2.set_ylim(ymin=0.,ymax=1.)
    plt.legend(lines,labels,ncol=1,loc='upper right') #, bbox_to_anchor=(1.0, 1.0),transform=ax1.transAxes)
#    fig.legend(lines,labels,ncol=1,loc='lower left', bbox_to_anchor=(1.0, 1.0))
    print 'save',river[point],re.split("_",pname[point])[2]+"_"+re.split("_",pname[point])[3],ylist[point],xlist[point],mean_obs[ylist[point],xlist[point]],std_obs[ylist[point],xlist[point]]
    #plt.savefig(assim_out+"/figures/sfcelv/"+river[point]+"_"+re.split("_",pname[point])[2]+"_"+re.split("_",pname[point])[3]+".png",dpi=300)
    plt.savefig(assim_out+"/figures/sfcelv/"+pname[point][2::]+".png",dpi=300)
    return 0
#
para_flag=1
#para_flag=0
#--
if para_flag==1:
    p=Pool(20)
    p.map(make_fig,np.arange(pnum))
    p.terminate()
else:
    map(make_fig,np.arange(pnum))