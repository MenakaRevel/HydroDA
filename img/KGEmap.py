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

#sys.path.append('../assim_out/')
import read_grdc as grdc
import cal_stat as stat
#import plot_colors as pc
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

#argvs = sys.argv

#experiment="E2O_HydroWeb22"
# experiment="VIC_BC_HydroWeb11"
# experiment="test_wse"
# experiment="DIR_WSE_E2O_HWEB_001"
# experiment="DIR_WSE_E2O_HWEB_002"
# experiment="ANO_WSE_E2O_HWEB_001"
# experiment="ANO_WSE_E2O_HWEB_002"
# experiment="ANO_WSE_E2O_HWEB_003"
# experiment="NOM_WSE_E2O_HWEB_001"
# experiment="NOM_WSE_E2O_HWEB_002"
# experiment="NOM_WSE_E2O_HWEB_003"
# experiment="NOM_WSE_E2O_HWEB_004"
# experiment="NOM_WSE_E2O_HWEB_005"
# experiment="NOM_WSE_E2O_HWEB_006"
# experiment="NOM_WSE_E2O_HWEB_007"
# experiment="NOM_WSE_E2O_HWEB_008"
# experiment="NOM_WSE_E2O_HWEB_009"
experiment="NOM_WSE_E2O_HWEB_010"

#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
# assim_out=pm.DA_dir()+"/out/"+experiment
assim_out="../out/"+experiment
print (assim_out)

# os.system("ln -sf "+assim_out+"/params.py params.py")
sys.path.append(assim_out)
import params as pm
#----
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    whereever the observed data contains nan
    
    this is used by all other functions, otherwise they will produce nan as 
    output
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    #mask = ~np.isnan(s) & ~np.isnan(o)
    #o_nonan = o[mask]
    #s_nonan = s[mask]

    #return o_nonan,s_nonan
    return data[:,0],data[:,1]
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
    # s = s[warmup+1:]
    # o = o[warmup+1:]
    # s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    B = np.mean(s) / np.mean(o)
    y = (np.std(s) / np.mean(s)) / (np.std(o) / np.mean(o))
    r = np.corrcoef(o, s)[0,1]
    return 1 - np.sqrt((r - 1) ** 2 + (B - 1) ** 2 + (y - 1) ** 2)
#----
def vec_par(LEVEL,ax=None):
    ax=ax or plt.gca()
    txt="KGEtmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec KGEtmp1.txt 1 "+str(LEVEL)+" > "+txt)
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
#----
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#----
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
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)
#----
mk_dir(assim_out+"/figures")
mk_dir(assim_out+"/figures/KGE")
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
rivermap=rivermap*(uparea>1e11)*1.0
#----
syear,smonth,sdate=2003,1,1 #spm.starttime()#2004#1991  2004,1,1 #
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
green ="green"
#-------
staid=[]
pname=[]
xlist=[]
ylist=[]
river=[]
rivernames = grdc.grdc_river_name_v396()
for rivername in rivernames:
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
    #print (rivername, grdc_id,station_loc)
    river.append([rivername]*len(station_loc))
    staid.append(grdc_id)
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)
#-----------------------
river=([flatten for inner in river for flatten in inner])
staid=([flatten for inner in staid for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
print (len(pname), len(xlist))
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
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
    fname=assim_out+"/assim_out/outflw/open/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=assim_out+"/assim_out/outflw/assim/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #-------------
    for point in np.arange(pnum):
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point])
        if ix2 == -9999 or iy2 == -9999:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]
        else:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]+opnfile[iy2,ix2]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]+asmfile[iy2,ix2]
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
cmap=mbar.colormap("H02")
cmap.set_under("w",alpha=0)
cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=0.0
vmax=1.0
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
os.system("./bin/txt_vector "+box+" "+pm.CaMa_dir()+" "+pm.mapname()+" > KGEtmp1.txt") 
#map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(2,10+1,1))
#--
for point in np.arange(pnum):
    org=grdc.grdc_dis(staid[point],syear,eyear-1)
    org=np.array(org)
    if np.sum(ma.masked_where(org!=-99.9,org))==0.0:
        print ("no obs", np.sum(ma.masked_where(org!=-99.9,org)))
        continue
    KGEasm=KGE(np.mean(asm[:,:,point],axis=1),org)
    KGEopn=KGE(np.mean(opn[:,:,point],axis=1),org)
    if KGEopn==1.00:
        print (KGEopn, staid[point], pname[point])
        continue
    KGEEAI=KGEasm #(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20) 
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
    c=cmapL(norm(KGEEAI))
    #print (lon,lat,NSEAI) #,NSEasm,NSEopn)
    # if NSEAI > 0.0:
    #     print lon,lat, "%3.2f %3.2f %3.2f"%(NSEAI, NSEasm, NSEopn)
    if KGEEAI >= 1e-20:
        ax.scatter(lon,lat,s=10,marker="o",edgecolors=c, facecolors=c,zorder=106)
    if KGEEAI < 1e-20:
        print staid[point], pname[point]
        #ax.scatter(lon,lat,s=10,marker="o",edgecolors="k", facecolors="k",zorder=106)
#--
cbar=m.colorbar(im,"right",size="2%",ticks=np.arange(vmin,vmax+0.001,0.2))
#plt.title(stitle)
plt.savefig(assim_out+"/figures/KGE/KGEscatter.png",dpi=300,bbox_inches="tight", pad_inches=0.05)
os.system("rm -r KGEtmp*.txt")