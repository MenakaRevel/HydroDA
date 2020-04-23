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
import my_colorbar as mbar

#assim_out="assim_out_E2O_womc"
#assim_out="assim_out_E2O_wmc"
#assim_out="assim_out_ECMWF_womc_baised_0"
#assim_out="assim_out_ECMWF_womc_baised"
#assim_out="assim_out_ECMWF_womc_baised_if"
#assim_out="assim_out_ECMWF_womc_baised_if_fixed1.10"
#assim_out="assim_out_ECMWF_womc_baised_if_adaptive"
#assim_out="assim_out"
#assim_out="assim_out_biased_womc"
##sys.path.append('../' +assim_out+'/')
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
mk_dir(assim_out+"/fig/NSE")
#----

#argvs = sys.argv

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

ratio=np.zeros((spix-npix)*(epix-wpix)).reshape([spix-npix,epix-wpix])
count=np.zeros((spix-npix)*(epix-wpix)).reshape([spix-npix,epix-wpix])
a2ones=np.ones((spix-npix)*(epix-wpix)).reshape([spix-npix,epix-wpix])

org=np.zeros([N,720,1440],np.float32)
asm=np.zeros([N,720,1440],np.float32)
#--
#fname = pm.CaMa_dir()+"/map/global_15min/rivout.bin"
#fname = pm.CaMa_dir()+"/map/glb_15min/outclm.bin"
fname = "../dat/mean_rivout_1960-2013.bin"
#fname="subriver.bin"
trueforo = np.fromfile(fname,np.float32).reshape([720,1440])
#trueforo = np.fromfile(fname,np.float32).reshape([2,720,1440])[0]
# ocean [ 0:ocean, 1:not ocean ]
#ocean = (trueforo[0,npix:spix,wpix:epix]<1e18) * 1
#----
rivnum="../dat/rivnum.bin"
rivnum = np.fromfile(rivnum,np.int32).reshape(720,1440)
#--
# river [ 0:not river, 1:river ]
river = (trueforo>100.)*(rivnum>0)*1.0
#river = (trueforo>0.)  * 1.0
#river = (trueforo[0,npix:spix,wpix:epix]>500.)  * 1

NSEs=["NSEopn","NSEasm","NSEAI"]
for NSE in NSEs:
    fname=assim_out+"/stat/"+NSE+".bin"
    NS=np.fromfile(fname,np.float32).reshape(720,1440)
    #--
    plt.close()
    #cmap=cm.viridis_r
    if NSE=="NSEAI":
        cmap=mbar.colormap("H01")
    else:
        cmap=mbar.colormap("H02")
    #cmap=cm.rainbow_r
    #cmap.set_under("w",alpha=0)
    NS=NS*river
    print NSE
    #print NS, np.shape(NS)

    resol=1
    plt.figure(figsize=(7*resol,3*resol))
    m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c')
    #m.drawcoastlines( linewidth=0.1, color='k' )
    m.fillcontinents(color=land,lake_color=water)
    m.drawmapboundary(fill_color=water)
    m.drawparallels(np.arange(south,north+0.1,20), labels = [1,0,0,0], fontsize=10,linewidth=0.1)
    m.drawmeridians(np.arange(west,east+0.1,40), labels = [0,0,0,1], fontsize=10,linewidth=0.1)
    #im = m.imshow(np.flipud(ratio),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
    data=NS[npix:spix,wpix:epix]
    #im = m.imshow(np.flipud(ma.masked_where(river[npix:spix,wpix:epix]!=1,data)),interpolation="nearest",vmin=1.0e-20, vmax=1,cmap=cmap,zorder=100)#
    im = m.imshow(np.flipud(ma.masked_less_equal(data,1.0e-20)),vmin=0.0, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
    cbar=m.colorbar(im,"right",size="2%")
    if NSE == "NSEAI":
        cbar.set_label("NSEAI")
        stitle="Nash-Sutcliff Efficency based Assimilation Index"
    else:
        cbar.set_label("NSE")
        stitle="Nash-Sutcliff Efficency"
    plt.title(stitle)
    plt.savefig(assim_out+"/fig/NSE/"+NSE+"_map.png",dpi=300,bbox_inches="tight", pad_inches=0.05)
