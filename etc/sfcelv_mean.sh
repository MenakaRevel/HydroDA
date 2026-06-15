#!/bin/bash
### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E40
#PBS -l select=1:ncpus=40:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N sfcelv_mean
#==
USER=`whoami`

# OMP Settings
NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/HydroDA/etc"

# expname
# expname="DIR_WSE_ECMWF_HWEB_011"
# expname="DIR_WSE_ECMWF_HWEB_012"
# expname="DIR_WSE_ECMWF_HWEB_013"
expname="DIR_WSE_ECMWF_HWEB_014"

# experiment directory
expdir="/cluster/data7/menaka/HydroDA/out/"$expname

# out dir
outdir="/cluster/data6/menaka/HydroDA/dat"

# out tag
tag="bias_corrupt_ECMWF_amz_06min_2009-2009"
# tag="bias_ECMWF_amz_06min_2009-2009"
# tag="corrupt_ECMWF_amz_06min_2009-2009"
# tag="49_ECMWF_amz_06min_2009-2009"

#CaMA-Flood directory
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# Map
mapname="amz_06min"

# syear, smon, sday
syear=2009
smon=1
sday=1

# eyear, emon, eday
eyear=2009
emon=12
eday=31

ens_num=49

N=`python calc_days.py $syear $smon $sday $eyear $emon $eday`

enum=1
while [ $enum -le $ens_num ];
do
    echo ./mean $syear $eyear $CaMa_dir $mapname $expdir $enum $outdir $tag $N
    time ./mean $syear $eyear $CaMa_dir $mapname $expdir $enum $outdir $tag $N &
    ## for parallel computation using multiple CPUs 
    NUM=`ps -U $USER | grep ./mean | wc -l | awk '{print $1}'`
    while [ $NUM -gt $NCPUS ];
    do
        sleep 1
        NUM=`ps -U $USER | grep ./mean | wc -l | awk '{print $1}'`
    done
    enum=$(( $enum + 1 ))
done

wait