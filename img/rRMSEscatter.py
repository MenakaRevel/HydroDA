#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import sys
import os
import calendar
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
import my_colorbar as mbar
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import os

#import plot_colors as pc
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

#argvs = sys.argv

#experiment="E2O_HydroWeb22"
# experiment="VIC_BC_HydroWeb11"
# experiment="test_wse"
# experiment="DIR_WSE_E2O_HWEB_001"
# experiment="DIR_WSE_E2O_HWEB_001"
# experiment="ANO_WSE_E2O_HWEB_001"
# experiment="ANO_WSE_E2O_HWEB_002"
# experiment="NOM_WSE_E2O_HWEB_001"
# experiment="NOM_WSE_E2O_HWEB_002"
# experiment="NOM_WSE_E2O_HWEB_003"
experiment="NOM_WSE_E2O_HWEB_004"
# experiment="NOM_WSE_E2O_HWEB_005"
# experiment="NOM_WSE_E2O_HWEB_006"
# experiment="NOM_WSE_E2O_HWEB_007"

#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
# assim_out=pm.DA_dir()+"/out/"+experiment
assim_out="../out/"+experiment
print (assim_out)

sys.path.append(assim_out)
# os.system("ln -sf ../gosh/params.py params.py")
import params as pm
import read_grdc as grdc
import read_hydroweb as hweb
import cal_stat as stat

conflag=pm.conflag()
print (conflag)
#====================================================================
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    where ever the observed data contains nan
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    return data[:,0],data[:,1]
#====================================================================
def vec_par(LEVEL,ax=None):
    ax=ax or plt.gca()
    txt="RMSEtmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec RMSEtmp1.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open(txt,"r")
    lines = f.readlines()
    f.close()
    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        # ix = int((lon1 + 180.)*(1/gsize))
        # iy = int((-lat1 + 90.)*(1/gsize))

        # if rivermap[iy,ix] == 0:
        #     continue

        if lon1-lon2 > 180.0:
            print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            print (lon1,lon2)
            lon2=-180.0
        #--------
        colorVal="w" 
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#====================================================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#====================================================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#====================================================================
def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    #s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s) 
    s,o = filter_nan(s,o)
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)
#==========================================================
def RMSE(s,o):
    """
    Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        RMSE: Root Mean Squre Error
    """
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    return np.sqrt(np.mean((s-o)**2))
#==========================================================
mk_dir(assim_out+"/figures")
mk_dir(assim_out+"/figures/rRMSE")
#----
fname=pm.CaMa_dir()+"/map/"+pm.mapname()+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
lon0   = float(filter(None, re.split(" ",lines[4]))[0])
lat0   = float(filter(None, re.split(" ",lines[7]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
#----
nextxy = pm.CaMa_dir()+"/map/"+pm.mapname()+"/nextxy.bin"
rivwth = pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivwth_gwdlr.bin"
rivhgt = pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivhgt.bin"
rivlen = pm.CaMa_dir()+"/map/"+pm.mapname()+"/rivlen.bin"
elevtn = pm.CaMa_dir()+"/map/"+pm.mapname()+"/elevtn.bin"
lonlat = pm.CaMa_dir()+"/map/"+pm.mapname()+"/lonlat.bin"
uparea = pm.CaMa_dir()+"/map/"+pm.mapname()+"/uparea.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
# #higher resolution data
# catmxy = pm.CaMa_dir()+"/map/"+pm.mapname()+"/1min/1min.catmxy.bin"
# catmxy = np.fromfile(catmxy,np.int16).reshape(2,ny*60,nx*60)
#----
rivnum="../dat/rivnum_"+pm.mapname()+".bin"
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#----
syear,smonth,sdate=pm.starttime() #2003,1,1 #2004#1991  2004,1,1 #
eyear,emonth,edate=pm.endtime() #2005,1,1 #2005,1,1 #
#month=1
#date=1
start_dt=datetime.date(syear,smonth,sdate)
end_dt=datetime.date(eyear,emonth,edate)
size=60

start=0
last=(end_dt-start_dt).days
N=int(last)

green2="greenyellow"
green ="green"
#-------
#-------
# mean & std for anomaly
mean_obss=np.zeros([pm.ens_mem(),ny,nx])
std_obss=np.zeros([pm.ens_mem(),ny,nx])
for num in np.arange(1,int(pm.ens_mem())+1):
    numch='%03d'%num
    fname=assim_out+"/assim_out/mean_sfcelv/meansfcelvC"+numch+".bin"
    mean_corr=np.fromfile(fname,np.float32).reshape([ny,nx])
    mean_obss[num-1,:,:]=mean_corr
    fname=assim_out+"/assim_out/mean_sfcelv/stdsfcelvC"+numch+".bin"
    std_corr=np.fromfile(fname,np.float32).reshape([ny,nx])
    std_obss[num-1,:,:]=std_corr
#-------
pname=[]
xlist=[]
ylist=[]
river=[]
EGM08=[]
EGM96=[]
rivernames  = ["AMAZONAS"]
for rivername in rivernames:
#   path = assim_out+"/figures/sfcelv/%s"%(rivername)
  station_loc,x_list,y_list,egm08,egm96 =hweb.get_hydroweb_loc(rivername,pm.mapname())
  river.append([rivername]*len(station_loc))
  pname.append(station_loc)
  xlist.append(x_list)
  ylist.append(y_list)
  EGM08.append(egm08)
  EGM96.append(egm96)
#-----------------------
river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
EGM08=([flatten for inner in EGM08 for flatten in inner])
EGM96=([flatten for inner in EGM96 for flatten in inner])
pnum=len(pname)
#########################################################
pnum=len(pname)
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
        #print (yyyy,mm,dd,numch)
#----------------------
# function to read data
def read_data(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    numch= inputlist[3]
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
    #fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=assim_out+"/assim_out/ens_xa/assim/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
    #fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #-------------
    for point in np.arange(pnum):
        ix1=xlist[point]
        iy1=ylist[point]
        tmp_opn[dt,num,point]=opnfile[iy1,ix1]
        tmp_asm[dt,num,point]=asmfile[iy1,ix1]        
        
#--------
p   = Pool(10)
res = p.map(read_data, inputlist)
opn = np.ctypeslib.as_array(shared_array_opn)
asm = np.ctypeslib.as_array(shared_array_asm)
p.terminate()		
#--------------------------------------------
land="#C0C0C0"
water="#FFFFFF"

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#cmap=make_colormap(colors_list)
# cmap=mbar.colormap("H01")
cmap=cm.coolwarm_r
# cmap.set_under("w",alpha=0)
cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=-0.2
vmax=0.2
norm=Normalize(vmin=vmin,vmax=vmax)
#------
# river width
sup=2
w=0.02
alpha=1
width=0.5
#------
print (west,east,south,north)
resol=1
fig=plt.figure()
G = gridspec.GridSpec(1,1)
ax=fig.add_subplot(G[0,0])
m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax=ax)
#m.drawcoastlines( linewidth=0.1, color='k' )
m.fillcontinents(color=land,lake_color=water,zorder=99)
#m.drawmapboundary(fill_color=water,zorder=100)
im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
im.set_visible(False)
m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
#--
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+pm.CaMa_dir()+" "+pm.mapname()+" > RMSEtmp1.txt") 
#map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(5,10+1,1))
#--
for point in np.arange(pnum):
    # time,org=hweb.HydroWeb_WSE(pname[point],syear,eyear-1)
    org=hweb.HydroWeb_continous_WSE(pname[point],syear,smonth,sdate,eyear-1,12,31,EGM08[point],EGM96[point])
    # print (np.sum((org!=-9999.0)*1.0))
    # if np.sum(ma.masked_where(org!=-9999.0,org))<0.0:
    if np.sum((org!=-9999.0)*1.0) == 0.0:
        # print ("no obs", np.sum((org!=-9999.0)*1.0))
        continue
    # org=np.array(org)+np.array(EGM08[point])-np.array(EGM96[point])
    # print (pname[point])
    # print ("=============================="+pname[point]+"====================================")
    #====================================
    org_mean=np.mean(ma.masked_equal(org,-9999.0))
    org_std=np.std(ma.masked_equal(org,-9999.0))
    if conflag==1:
        data0=ma.masked_where(org==-9999.0,org).filled(-9999.0)
    elif conflag==2:
        data0=ma.masked_where(org==-9999.0,(org-org_mean)).filled(-9999.0) #(org-mean_obs[ylist[point],xlist[point]])+mean_sfcelv[ylist[point],xlist[point]]
    elif conflag==3:
        data0=ma.masked_where(org==-9999.0,((org-org_mean)/(org_std+1e-20))).filled(-9999.0) #((org-mean_obs[ylist[point],xlist[point]])/(std_obs[ylist[point],xlist[point]]+1.0e-20))*std_sfcelv[ylist[point],xlist[point]]+mean_sfcelv[ylist[point],xlist[point]]
    #====================================
    opn0=np.zeros([N,pm.ens_mem()],np.float32)
    asm0=np.zeros([N,pm.ens_mem()],np.float32)
    for num in np.arange(0,int(pm.ens_mem())):
        if conflag==1:
            opn0[:,num]=opn[:,num,point]
            asm0[:,num]=asm[:,num,point]
        elif conflag==2:
            opn0[:,num]=opn[:,num,point]-mean_obss[num,ylist[point],xlist[point]]
            asm0[:,num]=asm[:,num,point]-mean_obss[num,ylist[point],xlist[point]]
        elif conflag==3:
            opn0[:,num]=(opn[:,num,point]-mean_obss[num,ylist[point],xlist[point]])/(std_obss[num,ylist[point],xlist[point]]+1e-20)
            asm0[:,num]=(asm[:,num,point]-mean_obss[num,ylist[point],xlist[point]])/(std_obss[num,ylist[point],xlist[point]]+1e-20)
    #====================================
    # print (np.shape(opn0))
    opn1=np.mean(opn0,axis=1)
    asm1=np.mean(asm0,axis=1)
    #====================================
    # print ("data0: ", data0, ma.masked_equal(data0,-9999.0).filled(0.0), np.mean(ma.masked_equal(data0,-9999.0)), conflag )
    RMSEasm=RMSE(asm1,data0)
    RMSEopn=RMSE(opn1,data0)
    #====================================
    if RMSEopn==0.0:
        print (RMSEopn, pname[point])
        continue
    rRMSE=(RMSEasm-RMSEopn)/(RMSEopn+1.0e-20)
    print (pname[point] , rRMSE, RMSEopn, RMSEasm) #, np.mean(asm1), np.mean(data0), 
    # np.sqrt(np.mean(ma.masked_where(data0==-9999.0,(data0-asm1)**2))) )
    # if rRMSE < 0.0:
    #     print (pname[point], rRMSE, RMSEopn, RMSEasm)
    #NSEAI=NSEasm
    ix=xlist[point]
    iy=ylist[point]
    if rivermap[iy,ix] !=1.0:
        continue
    #--------------
    #lon=lon0+ix*gsize
    #lat=lat0-iy*gsize
    lon=lonlat[0,iy,ix]
    lat=lonlat[1,iy,ix]
    c=cmapL(norm(rRMSE))
    #print (lon,lat,NSEAI) #,NSEasm,NSEopn)
    # if NSEAI > 0.0:
    #     print lon,lat, "%3.2f %3.2f %3.2f"%(NSEAI, NSEasm, NSEopn)
    ax.scatter(lon,lat,s=10,marker="o",edgecolors=c, facecolors=c,zorder=106)
    # if rRMSE >= 0.00:
    #     ax.scatter(lon,lat,s=10,marker="o",edgecolors=c, facecolors=c,zorder=106)
    # if rRMSE < 0.00:
    #     print staid[point], pname[point]
    #     ax.scatter(lon,lat,s=10,marker="d",edgecolors="k", facecolors="k",zorder=106)
#--
cbar=m.colorbar(im,"right",size="2%",ticks=np.arange(vmin,vmax+0.001,0.1),extend="both")
#plt.title(stitle)
plt.savefig(assim_out+"/figures/rRMSE/rRMSEscatter.png", dpi=300, bbox_inches="tight", pad_inches=0.05)
os.system("rm -r RMSEtmp*.txt")