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
from numpy import ma


#sys.path.append('../assim_out/')
os.system("ln -sf ../gosh/params.py params.py")
import params as pm
import read_grdc as grdc
import cal_stat as stat
#import plot_colors as pc
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

#argvs = sys.argv

experiment="E2O_HydroWeb4"
#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
assim_out=pm.DA_dir()+"/out/"+experiment
print assim_out
#assim_out="assim_out_E2O_wmc"
#assim_out="assim_out_E2O_womc_0"
#assim_out="assim_out_ECMWF_womc_baised_0"
#assim_out="assim_out_ECMWF_womc_baised"
#assim_out="assim_out_ECMWF_womc_baised_if"
#assim_out="assim_out_ECMWF_womc_baised_0.90"
#assim_out="assim_out_ECMWF_womc_baised_0.80"
#assim_out="assim_out_ECMWF_womc_baised_0.75"
#assim_out="assim_out_ECMWF_womc_baised_0.20"
#assim_out="assim_out_ECMWF_womc_baised_0.50"
#assim_out="assim_out_ECMWF_womc_baised_if_fixed1.08"
#assim_out="assim_out_ECMWF_womc_baised_if_fixed1.10"
#assim_out="assim_out_ECMWF_womc_baised_if_adaptive"
#assim_out="assim_out"
#assim_out="assim_out_biased_womc"
#assim_out="assim_out_biased_wmc"
#os.system("mkdir ../assim_out/img")
#os.system("mkdir ../assim_out/img/disgraph")
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
mk_dir(assim_out+"/figures/disgraph")
#----
year=2004#1991
month=1
date=1
start_dt=datetime.date(year,month,date)
size=60

start=0
last=365#int(argvs[1])
if calendar.isleap(year):
    last=366
else:
    last=365

#last=89
N=int(last)

green2="greenyellow"
green ="green"
#colors = pc.color_assimd()
staid=[]
pname=[]
xlist=[]
ylist=[]
river=[]
#--
rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZON","MEKONG","IRRAWADDY","VOLGA", "NIGER","YUKON","DANUBE"] #,"INDUS"] #["AMAZONAS"]#["CONGO"]#
#rivernames = grdc.grdc_river_name_v396()
#rivernames = grdc.grdc_river_name()
#for rivername in rivernames:
#  #station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
#  station_loc,x_list,y_list = grdc.get_grdc_station_v396(rivername)
#  #print rivername, station_loc
#  for station in station_loc:
#    gid=grdc.get_id(station)
#    if gid== -9999:
#      continue
#    path = "../"+assim_out+"/figures/disgraph/%s"%(rivername)
#    #print path
#    mk_dir(path)
#    ix, iy = grdc.get_loc_v394(gid)
#    print station,gid, ix ,iy
#    river.append(rivername)
#    pname.append(station)
#    xlist.append(ix)
#    ylist.append(iy)
for rivername in rivernames:
  path = assim_out+"/figures/disgraph/%s"%(rivername)
  print path
  mk_dir(path)
  #station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
  grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
  print rivername, grdc_id,station_loc
  river.append([rivername]*len(station_loc))
  staid.append(grdc_id)
  pname.append(station_loc)
  xlist.append(x_list)
  ylist.append(y_list)
#--
#for rivername in ["LENA","NIGER","INDUS","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS"]:
#for rivername in ["MEKONG","AMAZONAS"]:
#  if rivername=="LENA":
#      pname.append(["L1","L2","L3"])
#      xlist.append([1233,1218,1206])
#      ylist.append([  71,  98, 117])
#      river.append([rivername,rivername,rivername])
#  if rivername=="NIGER":
#      pname.append(["N1","N2","N3","N4","N5","N6"])
#     xlist.append([744,744,732,712,704,700])
#      ylist.append([342,324,310,292,295,303])
#     river.append([rivername,rivername,rivername,rivername,rivername,rivername])
# if rivername=="AMAZONAS":
#      pname.append(["B","E","F","G"])
#      xlist.append([515,447,464,420])
#      ylist.append([364,366,416,367])
#      river.append([rivername,rivername,rivername,rivername])
#  if rivername=="MEKONG":
#      pname.append(["Me_D","Me_M","Me_U"])
#      xlist.append([1143,1139,1127])
#      ylist.append([ 319, 293, 282])
#      river.append([rivername,rivername,rivername])
#  if rivername=="MISSISSIPPI":
#      pname.append(["Mi_D","Mi_M","Mi_U","Mi_U2"]) #最後はミズーリ川
#      xlist.append([361,362,345,308])
#      ylist.append([241,214,137,167])
#      river.append([rivername,rivername,rivername,rivername])
#  if rivername=="OB":
#      pname.append(["O_D","O_M","O_U"])
#      xlist.append([995,996,1048])
#      ylist.append([ 92,121, 159])
#      river.append([rivername,rivername,rivername])
#  if rivername=="CONGO":
#      pname.append(["C_D","C_M","C_U"])#,"C1","C2","C3","C4","C5","C6","C7","C8"])
#      xlist.append([772,813,834])#,769,784,789,809,821,824,802,794])
#      ylist.append([383,353,397])#,383,372,363,343,359,376,379,342])
#      river.append([rivername,rivername,rivername])#,rivername,rivername,rivername,rivername,rivername,rivername,rivername,rivername])
#  if rivername=="INDUS":
#      pname.append(["I_D","I_Sub","I_M"])
#      xlist.append([992,995,1003])
#      ylist.append([251,252,233])
#      river.append([rivername,rivername,rivername])
#  if rivername=="CONGO":
#    river.append([rivername]*12)
#    pname.append(["C4","C1","C10","C9","C7","C5","C3","C2","C8","C6","C11","C12"])
#    xlist.append([813,834,772,784,789,809,821,824,802,794,808,806])
#    ylist.append([353,397,383,372,363,343,359,376,379,342,351,351])

river=([flatten for inner in river for flatten in inner])
staid=([flatten for inner in staid for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
print len(pname), len(xlist)
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])


pnum=len(pname)

org=[]
opn=[]
asm=[]

swt={}
for point in np.arange(pnum):
    swt[point] = []

for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    print yyyy,mm,dd

#    fname="../data/mesh_day%02d.bin"%(SWOT_day(yyyy,mm,dd))
#    mesh_in=np.fromfile(fname,np.float32).reshape([640,1440])
#    mesh=(mesh_in>=10)*(mesh_in<=60)
#    meshP=mesh-1000*(mesh<0.1)

#    fname="../sat/observation_day%02d.bin"%(SWOT_day(yyyy,mm,dd))
#    meshP=np.fromfile(fname,np.int32).reshape([720,1440])
#    #meshP=(meshP>=1)*1


#    # make org
#    fname=assim_out+"/rivout/true/rivout"+yyyy+mm+dd+".bin"
#    orgfile=np.fromfile(fname,np.float32).reshape([720,1440])
#
#    org_frag=[]
#    for point in np.arange(pnum):
#        #print point
#        xpoint=xlist[point]
#        ypoint=ylist[point]
#        org_frag.append(orgfile[ypoint,xpoint])
#        #---SWOT--
#        #if meshP[ypoint-40,xpoint] >= 1:
#        if meshP[ypoint,xpoint] >= 1:
#          if point not in swt.keys():
#            swt[point] = [day]
#          else:
#            swt[point].append(day)
#
#    org.append(org_frag)

    # make asm and opn
    opn_ens=[]
    asm_ens=[]
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num

        fname=assim_out+"/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
        #fname="../CaMa_out/"+yyyy+mm+dd+"C"+numch+"/rivout"+yyyy+".bin"
        opnfile=np.fromfile(fname,np.float32).reshape([720,1440])

        fname=assim_out+"/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
        #fname="../CaMa_out/"+yyyy+mm+dd+"A"+numch+"/rivout"+yyyy+".bin"
        asmfile=np.fromfile(fname,np.float32).reshape([720,1440])

        opn_frag=[]
        asm_frag=[]
        for point in np.arange(pnum):
            xpoint=xlist[point]
            ypoint=ylist[point]
            opn_frag.append(opnfile[ypoint,xpoint])
            asm_frag.append(asmfile[ypoint,xpoint])

        opn_ens.append(opn_frag)
        asm_ens.append(asm_frag)


    opn.append(opn_ens)
    asm.append(asm_ens)

#org=np.array(org,dtype=np.float32)
opn=np.array(opn,dtype=np.float32)
asm=np.array(asm,dtype=np.float32)

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
    fig, ax1 = plt.subplots()
    org=grdc.grdc_dis(staid[point],year,year)
    org=np.array(org)
    ax1.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="black",linewidth=0.7,zorder=101) #,marker = "o",markevery=swt[point])
#    ax1.plot(np.arange(start,last),hgt[:,point],label="true",color="gray",linewidth=0.7,linestyle="--",zorder=101)
#    plt.plot(np.arange(start,last),org[:,point],label="true",color="black",linewidth=0.7)

    for num in np.arange(0,pm.ens_mem()):
        ax1.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color="blue",linewidth=0.1,alpha=0.3,zorder=102)
        ax1.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color="red",linewidth=0.1,alpha=0.3,zorder=103)
#        plt.plot(np.arange(start,last),opn[:,num,point],label="corrupted",color="blue",linewidth=0.3,alpha=0.5)
#        plt.plot(np.arange(start,last),asm[:,num,point],label="assimilated",color="red",linewidth=0.3,alpha=0.5)
    # draw mean of ensembles
    ax1.plot(np.arange(start,last),np.mean(opn[:,:,point],axis=1),label="corrupted",color="blue",linewidth=0.8,alpha=1,zorder=104)
    ax1.plot(np.arange(start,last),np.mean(asm[:,:,point],axis=1),label="assimilated",color="red",linewidth=0.8,alpha=1,zorder=106)

#    plt.ylim(ymin=)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('discharge (m$^3$/s)', color='k')
    ax1.set_xlim(xmin=0,xmax=last+1)
    ax1.tick_params('y', colors='k')
    # scentific notaion
    ax1.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    xxlist=np.linspace(15,N-15,12)
    xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    ax1.set_xticks(xxlist)
    ax1.set_xticklabels(xxlab,fontsize=10)

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
    #ax2.plot(np.arange(start,last),ma.masked_less_equal(aivn,1.0e-20),color="green",zorder=104,marker = "o", alpha=0.5,linewidth=0.5,markevery=swt[point])
#    ax2.plot(np.arange(start,last),ma.masked_where(error<=0.1,aivn),color=green,marker = "o",markeredgecolor =green, alpha=1,linewidth=0.5,markevery=swt[point],markersize=3,zorder=100)#,label="AI")
#    ax2.plot(np.arange(start,last),ma.masked_where(error>0.1,aivn),color=green2,marker = "o",markeredgecolor =green2, alpha=1,linewidth=0.5,markevery=swt[point],markersize=3,zorder=100)#,label="AI")
#    ax2.set_ylabel('AI', color='green')
#    ax2.tick_params('y', colors='green')
#    ax2.set_ylim(ymin=0.,ymax=1.)
#    ax2.set_xlim(xmin=0,xmax=last+1)
#    print swt[point]
    print 'save',river[point],pname[point]
    plt.savefig(assim_out+"/figures/disgraph/"+river[point]+"/"+pname[point]+".png",dpi=500)
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
#para_flag=0
#--
if para_flag==1:
    p=Pool(20)
    p.map(make_fig,np.arange(pnum))
    p.terminate()
else:
    map(make_fig,np.arange(pnum))
