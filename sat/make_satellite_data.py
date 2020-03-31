#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import sys
import math
import itertools
from multiprocessing import Pool
from multiprocessing import Process
############################################
def calc_dis(latlonlist):
    lat=latlonlist[0]
    lon=latlonlist[1]
    diff=np.sqrt((lat-orbit[:,1])**2+((lon-orbit[:,2])*math.cos(lat*3.1416/180.))**2)
    #print diff

    min_diff=min(diff[orbit[:,0]==day])
    #print diff,min_diff
    min_index=int(min(orbit[:,3][diff==min_diff]))
    #print min_diff-diff[min_index]

    #距離を測る
    # 最初と最後の場合の例外処理の追加が必要
    # Measure Distance
    # Need to add exception handaling for first and last cases
    if min_index==61068:
        min_index=min_index-1
    if min_index==0:
        min_index=min_index+1

    x0=lon
    y0=lat
    x1=orbit[min_index,2]
    y1=orbit[min_index,1]
    x2=orbit[min_index+1,2]
    y2=orbit[min_index+1,1]
    if x2==x1 and y2==y1:
        x2=orbit[min_index+2,2]
        y2=orbit[min_index+2,1]
    x3=orbit[min_index-1,2]
    y3=orbit[min_index-1,1]
    if x3==x1 and y3==y1:
        x3=orbit[min_index-2,2]
        y3=orbit[min_index-2,1]

    a1=(y2-y1)/((x2-x1+1e-20)*math.cos(lat*3.1416/180.))
    a2=(y3-y1)/((x3-x1+1e-20)*math.cos(lat*3.1416/180.))

    a=(a1+a2)/2.
    b=-1.
    c=0.  #-(a1+a2)/2.*x1+y1
    # すぐ上のはx1,y1を座標軸としたので0でOKになった

    d=abs(a*(x0-x1)*math.cos(lat*3.1416/180.)+b*(y0-y1)+c)/math.sqrt(a**2+b**2)
    #print d,(d>=10 and d<=60)

    #緯度方向、経度方向の距離に
    ds=-2./(a1+a2)
    dsA=d*ds/math.sqrt(1+ds**2) #緯度方向
    A=d/math.sqrt(1+ds**2) #経度方向

    #kmに直す
    eqa_len=40076.954 #km
    mer_len=40009.000 #km

    mer_dis=dsA*mer_len/360. #緯度方向
    eqa_dis=A*eqa_len/360. #経度方向

    dis=math.sqrt(mer_dis**2+eqa_dis**2) #distance
    if lat%10==0 and lon%10==0:
        print lat,lon#,dis
    return dis
###################################################################
def observation_error():
    """observation error of WSE depending on the L*W of each pixel
    used sigma*(1/l)*(1/w) l=k*L, w=q*W  Rodrigaz et al 2017:
    According to CaMa k=0.25, q=0.85"""
    k=1.00 # assume nearest part to the unit catchment
    q=1.00 # used 1.0 -> river width variability is 30%
    rivlen=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivlen.bin",np.float32).reshape(720,1440)
    rivwth=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivwth_gwdlr.bin",np.float32).reshape(720,1440)
    nextx=(np.fromfile(pm.CaMa_dir()+"/map/glb_15min/nextxy.bin",np.int32).reshape(2,720,1440)[0]!=-9999)*1.0
    rivlen=1.0 #rivlen*1.0e-3 #used as one kilmeter
    rivwth=rivwth*1.0e-3
    area=(k*rivlen)*(q*rivwth)
    obs_err=pm.ovs_err()*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx
    #obs_err=pm.ovs_err()*(1/(area+1.0e-20))*nextx
    # if water area < 1.0 km2 -> 0.25
    obs_err=obs_err*(area>=1.0)*1.0+0.25*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx*(area<1.0)*1.0
    obs_err=ma.masked_where(area<0.625,obs_err).filled(0.25) # 25cm for < 1km^2 water area
    obs_err=obs_err*nextx
    obs_err=obs_err.astype(np.float32)
    obs_err.tofile("./obs_err.bin")
    return 0
#####################################################################
# main code
f = open('./SWOT_Science_sept2015_Nadir.kml')
data=f.readlines()
f.close()

incoord=0 #0は外 1はin

orbit=[] #[day,lat,lon,index] 4 x index数 次元

index=0

for line in data:
    #dayの特定
    # Dayって文字を探して次の２ケタを拾う
    if line.find("Day")>=0:
        day=int(line[line.find("Day")+3:line.find("Day")+5])

    if line.find("</coordinates>")>=0:
        incoord=0

    #読み込み
    if incoord==1:
        lon=float(line[:line.find(",")])
        left=line[line.find(",")+1:]
        lat=float(left[:left.find(",")])
        orbit.append([day,lat,lon,index])
        index=index+1

    if line.find("<coordinates>")>=0:
        incoord=1

orbit=np.asarray(orbit)

# dayを入力
for day in np.arange(1,22,1):
    # あるinput lat,lonに対して
    print day
    latlist=np.arange(80.,-80.,-0.25)
    lonlist=np.arange(-180.,180.,0.25)

    latlonlist=list(itertools.product(latlist,lonlist))

    p = Pool()
    outmesh=p.map(calc_dis,latlonlist)

    day_p="%02d"%day
    wfile="./mesh_day"+day_p+".bin"

    outmesh=np.asarray(outmesh)
    outmesh=outmesh.astype(np.float32)
    outmesh.tofile(wfile)

    print outmesh

#
#
## メッシュを作る
#lon_mesh=np.asarray([[np.arange(-180.,180.,0.25)] for i in np.arange(180*4)])
#lat_mesh=np.asarray([[np.ones(360*4)*i] for i in np.arange(90.,-90.,-0.25)])
#
#diff=lon_mesh
#
#for i=np.arange(-90.,90.,0.25): #lat
#    for j=np.arange(-180.,180.,0.25): #lon
#        diff.append(sqrt((i-lat)**2+(j-lon)**2))
#
#min(diff)

#make observation error variance
observation_error()
