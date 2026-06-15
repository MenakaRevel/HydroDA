#!/opt/local/bin/pyuhon
#-*- coding: utf-8 -*-
 
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
from mpl_toolkits.basemap import Basemap
import os
import calendar


#assim_out="assim_out_E2O_womc"
#assim_out="assim_out_E2O_wmc"
#assim_out="assim_out_biased_womc"
#assim_out="assim_out_ECMWF_womc_baised_if_fixed1.10"
#sys.path.append('../'+assim_out+'/')
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
mk_dir(assim_out+"/fig/AI")
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
#--
# run calc_stat.py to create the statistical maps
ratio=np.fromfile(assim_out+"/stat/annualmeanAI.bin",np.float32).reshape(720,1440)
###ratio=np.zeros((spix-npix)*(epix-wpix)).reshape([spix-npix,epix-wpix])
###count=np.zeros((spix-npix)*(epix-wpix)).reshape([spix-npix,epix-wpix])
#*******************************************************
###if os.path.exists("../"+assim_out+"/img/AImap/annualmeanAI.bin"):
###    ratio=np.fromfile("../"+assim_out+"/img/AImap/annualmeanAI.bin",np.float32).reshape(720,1440)
###else:
###    for day in np.arange(0,lastday):
###        #for day in np.arange(100,110):
###        # analyse date
###        target_dt=start_dt+datetime.timedelta(days=day)
###        yyyy='%04d' % (target_dt.year)
###        mm='%02d' % (target_dt.month)
###        dd='%02d' % (target_dt.day)
###        print yyyy,mm,dd
###
###        # next day name
###        next_dt=start_dt+datetime.timedelta(days=day+1)
###        nxt_yyyy='%04d' % (next_dt.year)
###        nxt_mm='%02d' % (next_dt.month)
###        nxt_dd='%02d' % (next_dt.day)
###
###        # True Discharge
###        fname="../assim_out/rivout/true/rivout"+yyyy+mm+dd+".bin"
###        org=np.fromfile(fname,np.float32).reshape([720,1440])
###
###        # open loop
###        opn=[]
###        for num in np.arange(1,pm.ens_mem()+1):
###            numch = "%03d" % num
###            fname = "../assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
###            opn.append(np.fromfile(fname,np.float32).reshape([720,1440]))
###        opn = np.array(opn)
###        opn_mean=np.mean(opn,axis=0)
###
###        # assimilated
###        asm=[]
###        for num in np.arange(1,pm.ens_mem()+1):
###            numch = "%03d" % num
###            fname = "../assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
###            asm.append(np.fromfile(fname,np.float32).reshape([720,1440]))
###        asm = np.array(asm)
###        asm_mean=np.mean(asm,axis=0)
###
###        # assimilation index 計算
###        #ai=1-abs((asm_mean[npix:spix,wpix:epix]-opn_mean[npix:spix,wpix:epix])/(org[npix:spix,wpix:epix]-opn_mean[npix:spix,wpix:epix]+1e-20)-1)
###        ai=1.- np.absolute((asm_mean[npix:spix,wpix:epix]-opn_mean[npix:spix,wpix:epix])/((org[npix:spix,wpix:epix]-opn_mean[npix:spix,wpix:epix])+1e-20)-1.)
###    # read restart file for making ocean mask
###        #fname = "../CaMa_in/restart/true/restart" + nxt_yyyy + nxt_mm + nxt_dd + "T.bin"
###        fname = pm.CaMa_dir()+"/map/global_15min/rivout.bin"
###        trueforo = np.fromfile(fname,np.float32).reshape([2,720,1440])
###        # ocean [ 0:ocean, 1:not ocean ]
###        #ocean = (trueforo[0,npix:spix,wpix:epix]<1e18) * 1
###
###        # river [ 0:not river, 1:river ]
###        river = (trueforo[0,npix:spix,wpix:epix]>500.)  * 1
###
###        # error < 10%
###        error=((np.absolute(org[npix:spix,wpix:epix]-opn_mean[npix:spix,wpix:epix])/(org[npix:spix,wpix:epix]+1e-20))>0.1)*(1)
###        error=np.nan_to_num(error)
###        #--
###        river = river*error
###    
###        # ratio
###        ratio_n = ai * river #* ocean
###        ratio_n = np.ma.fix_invalid(ratio_n,fill_value=0.0)
###        ratio = ratio + (ratio_n<0)*0+(ratio_n>1)*1+(ratio_n>=0)*(ratio_n<=1)*ratio_n
###        count = count + river #* ocean
###    
###    ratio = ratio / (count.astype(np.float32)+1.0e-20)
###    river = (trueforo[0,npix:spix,wpix:epix]>500.)  * 1
###    ratio = ratio * river
###    
###    ai = ai * river
# assim/true
plt.close()
cmap=cm.viridis_r
cmap.set_under("w",alpha=0)

resol=1
plt.figure(figsize=(7*resol,3*resol))
m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c')
#m.drawcoastlines( linewidth=0.3, color='k' )
m.fillcontinents(color=land,lake_color=water)
m.drawmapboundary(fill_color=water)
m.drawparallels(np.arange(south,north+0.1,20), labels = [1,0,0,0], fontsize=10,linewidth=0.1)
m.drawmeridians(np.arange(west,east+0.1,40), labels = [0,0,0,1], fontsize=10,linewidth=0.1)
fname = pm.CaMa_dir()+"/map/glb_15min/outclm.bin"
trueforo = np.fromfile(fname,np.float32).reshape([2,720,1440])[0]
river=(trueforo>100.)*1.0
ratio=np.ma.fix_invalid(ratio).data
ratio=ratio*river
data=ratio[npix:spix,wpix:epix]
im = m.imshow(np.flipud(data),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
#im = m.imshow(np.flipud(ai),vmin=1e-20, vmax=1,interpolation="nearest",cmap=cmap,zorder=100)
cbar=m.colorbar(im,"right",size="2%")
cbar.set_label("annual mean AI")
plt.title("annual mean Assimilation Index ")#+yyyy+"-"+mm+"-"+dd)
plt.savefig(assim_out+"/fig/AI/AImap.png",dpi=300,bbox_inches="tight", pad_inches=0.05)

