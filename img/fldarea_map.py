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
# Link the params.py in the experiment dir
# os.system("ln -sf ../gosh/params_real.py params.py")
# import params as pm

#import plot_colors as pc
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

# os.system("rm -rf params.py")
#argvs = sys.argv

# experiment="E2O_HydroWeb23"
# experiment="VIC_BC_HydroWeb11"
# experiment="test_wse"
# experiment="test_virtual"
# experiment="DIR_WSE_E2O_HWEB_001"
# experiment="DIR_WSE_E2O_HWEB_002"
# experiment="DIR_WSE_E2O_HWEB_003"
# experiment="DIR_WSE_E2O_HWEB_004"
# experiment="ANO_WSE_E2O_HWEB_001"
# experiment="ANO_WSE_E2O_HWEB_002"
# experiment="ANO_WSE_E2O_HWEB_003"
# experiment="ANO_WSE_E2O_HWEB_004"
# experiment="NOM_WSE_E2O_HWEB_001"
# experiment="NOM_WSE_E2O_HWEB_002"
experiment="NOM_WSE_E2O_HWEB_003"
# experiment="NOM_WSE_E2O_HWEB_004"
# experiment="NOM_WSE_E2O_HWEB_005"
# experiment="NOM_WSE_E2O_HWEB_006"
# experiment="NOM_WSE_E2O_HWEB_007"
# experiment="NOM_WSE_E2O_HWEB_008"
# experiment="NOM_WSE_E2O_HWEB_009"
# experiment="NOM_WSE_E2O_HWEB_010"
# experiment="NOM_WSE_E2O_HWEB_011"
# experiment="NOM_WSE_E2O_HWEB_012"
# experiment="NOM_WSE_E2O_HWEB_013"

#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
# assim_out=pm.DA_dir()+"/out/"+experiment
# assim_out="../out/"+experiment
assim_out="/cluster/data7/menaka/HydroDA/out/"+experiment
print (assim_out)

sys.path.append(assim_out)
import params as pm
import read_grdc as grdc
import cal_stat as stat
#====================================================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#==========================================================
mk_dir(assim_out+"/figures")
mk_dir(assim_out+"/figures/flood")
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
syear,smonth,sdate=pm.starttime()#2004#1991  2004,1,1 # 2003,1,1 #
eyear,emonth,edate=2010,1,1 #pm.endtime() #2005,1,1 # 2005,1,1 # 2010,1,1 #
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
org=[]
opn=[]
asm=[]

# multiprocessing array
# result       = np.ctypeslib.as_ctypes(np.zeros((size, size)))
# shared_array = sharedctypes.RawArray(result._type_, result)
opn=np.ctypeslib.as_ctypes(np.zeros([N,pm.ens_mem(),ny,nx],np.float32))
shared_array_opn  = sharedctypes.RawArray(opn._type_, opn)
asm=np.ctypeslib.as_ctypes(np.zeros([N,pm.ens_mem(),ny,nx],np.float32))
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

def read_data(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    numch= inputlist[3]
    # print (yyyy,mm,dd,numch)
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
    fname=assim_out+"/assim_out/fldarea/open/fldarea"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=assim_out+"/assim_out/fldarea/assim/fldarea"+yyyy+mm+dd+"_"+numch+".bin"
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #-------------
    tmp_opn[dt,num,:,:]=opnfile #[iy1,ix1]
    tmp_asm[dt,num,:,:]=asmfile #[iy1,ix1]
    
#--------
p   = Pool(20)
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
# cmap=cm.get_cmap("bwr")
# cmap=cm.get_cmap("Spectral_r")
# cmap=cm.get_cmap("PiYG")
# cmap=cm.get_cmap("PRGn")
# cmap=mbar.colormap("H02")
cmap=cm.get_cmap("ocean_r")
# cmap=cm.get_cmap("YlGnBu")
# cmap=cm.get_cmap("nipy_spectral_r")
cmap.set_under("w",alpha=0)
cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=0.0
vmax=1.e2#max(np.amax(asm),np.amax(opn))
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
G = gridspec.GridSpec(1,3)
ax1=fig.add_subplot(G[0,0])
m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax=ax1)
# m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)
m.fillcontinents(color="gray",lake_color=water,zorder=99)
#m.drawmapboundary(fill_color=water,zorder=100)
data=ma.masked_greater(np.mean(asm,axis=(0,1)),1.0e10)*1.0e-6
data=ma.masked_less_equal(data,0.0)
im=m.imshow(data,origin="upper",interpolation="nearest",cmap=cmapL,vmin=vmin,vmax=vmax,zorder=110)#,norm=norm,zorder=101)
print (data)
m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=4,linewidth=0,zorder=102)
m.drawmeridians(np.arange(west,east+0.1,10), labels = [0,0,0,1], fontsize=4,linewidth=0,zorder=102)
#--
ax2=fig.add_subplot(G[0,1])
m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax=ax2)
# m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)
m.fillcontinents(color="gray",lake_color=water,zorder=99)
#m.drawmapboundary(fill_color=water,zorder=100)
data=ma.masked_greater(np.mean(opn,axis=(0,1)),1.0e10)*1.0e-6
data=ma.masked_less_equal(data,0.0)
im=m.imshow(data,origin="upper",interpolation="nearest",cmap=cmapL,vmin=vmin,vmax=vmax,zorder=110)#,norm=norm,zorder=101)

m.drawparallels(np.arange(south,north+0.1,5), labels = [0,0,0,0], fontsize=4,linewidth=0,zorder=102)
m.drawmeridians(np.arange(west,east+0.1,10), labels = [0,0,0,1], fontsize=4,linewidth=0,zorder=102)
#--
cbar=m.colorbar(im,"right",size="2%",ticks=np.arange(vmin,vmax+0.001,100),extend="max")
# cbar.set_label("$\Delta$$correlation$") #$(r_{assim} - r_{open})$")
# cbar.set_label("$r $ $correlation$") #$(r_{assim} - r_{open})$")
#plt.title(stitle)
plt.savefig(assim_out+"/figures/flood/fldarea.png",dpi=500,bbox_inches="tight", pad_inches=0.05)