#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
from mpl_toolkits.basemap import Basemap
import os
from numpy import ma
import calendar

#assim_out="assim_out_E2O_womc"
assim_out="assim_out_E2O_wmc"
#assim_out="assim_out_ECMWF_womc_baised_0"
#assim_out="assim_out_ECMWF_womc_baised_if_fixed1.10"
#assim_out="assim_out"
#assim_out="assim_out_biased_womc"
os.system("ln -sf ../gosh/params.py params.py")
import params as pm

experiment="E2O_wmc_06"
#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
print assim_out
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#----
mk_dir(assim_out+"/fig")
mk_dir(assim_out+"/fig/KGE")
#----

year=2004
month=1
date=1
start_dt=datetime.date(year,month,date)
size=60

south= -90
north=  90
west= -180
east=  180

land="#FFFFFF"
water="#C0C0C0"

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#lastday=int(argvs[1])
lastday=365 #int(argvs[1])
if calendar.isleap(year):
   lastday=366
else:
   lastday=365
N=lastday

#--
fname = "../dat/mean_rivout_1960-2013.bin"
trueforo = np.fromfile(fname,np.float32).reshape([720,1440])
#----
rivnum="../dat/rivnum.bin"
rivnum = np.fromfile(rivnum,np.int32).reshape(720,1440)
#--
# river [ 0:not river, 1:river ]
river = (trueforo>100.)*(rivnum>0)*(rivnum<=1000)*1.0

# run  calc_stat.py to create NSEasm.bin, NSEopn.bin, and NSEAI
#*******************************************************
KGEs=["KGEopn","KGEasm"]
#if os.path.exists("../"+assim_out+"/img/AImap/"+KGE+".bin"):
for KGE in KGEs:
    KG=np.fromfile(assim_out+"/stat/"+KGE+".bin",np.float32).reshape(720,1440)
    print KGE
    KG=KG*river
    #print KG, np.shape(KG)
    # assim/true
    plt.close()
    #cmap=cm.viridis_r
    cmap=cm.rainbow_r
    #cmap.set_under("w",alpha=0)
    resol=1
    plt.figure(figsize=(7*resol,3*resol))
    m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c')
    #m.drawcoastlines( linewidth=0.3, color='k' )
    m.fillcontinents(color=land,lake_color=water)
    m.drawmapboundary(fill_color=water)
    m.drawparallels(np.arange(south,north+0.1,20), labels = [1,0,0,0], fontsize=10,linewidth=0.1)
    m.drawmeridians(np.arange(west,east+0.1,40), labels = [0,0,0,1], fontsize=10,linewidth=0.1)
    #im = m.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
    data=KG[npix:spix,wpix:epix]
    im = m.imshow(np.flipud(ma.masked_where(river[npix:spix,wpix:epix]<=0,data)),vmin=0, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)#
    cbar=m.colorbar(im,"right",size="2%")
    cbar.set_label("KGE")
    plt.title("Kling Gupta Efficiency")#+yyyy+"-"+mm+"-"+dd)
    plt.savefig(assim_out+"/fig/KGE/"+KGE+"_map.png",dpi=300,bbox_inches="tight", pad_inches=0.05)


KGEasm=np.fromfile(assim_out+"/stat/KGEasm.bin",np.float32).reshape(720,1440)
KGEopn=np.fromfile(assim_out+"/stat/KGEopn.bin",np.float32).reshape(720,1440)

# KGE differnce 
plt.close()
#cmap=cm.viridis_r
cmap=cm.bwr_r
#cmap.set_under("w",alpha=0)

resol=1
plt.figure(figsize=(7*resol,3*resol))
m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c')
#m.drawcoastlines( linewidth=0.3, color='k' )
m.fillcontinents(color=land,lake_color=water)
m.drawmapboundary(fill_color=water)
m.drawparallels(np.arange(south,north+0.1,20), labels = [1,0,0,0], fontsize=10,linewidth=0.1)
m.drawmeridians(np.arange(west,east+0.1,40), labels = [0,0,0,1], fontsize=10,linewidth=0.1)
#im = m.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
data=KGEasm[npix:spix,wpix:epix]-KGEopn[npix:spix,wpix:epix]
im = m.imshow(np.flipud(ma.masked_where(river[npix:spix,wpix:epix]<=0,data)),vmin=-0.5, vmax=0.5,interpolation="nearest",cmap=cmap,zorder=100)#
cbar=m.colorbar(im,"right",size="2%",extend="both",ticks=np.arange(-0.5,0.5+0.1,0.25))
cbar.set_label("KGE differnce")
plt.title("Kling Gupta Efficiency")#+yyyy+"-"+mm+"-"+dd)
plt.savefig(assim_out+"/fig/KGE/KGEdiff_map.png",dpi=300,bbox_inches="tight", pad_inches=0.05)