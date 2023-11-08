#!/bin/bash
#***********************************************************************************************
# Data Assimilation using LETKF and empircal local patches [Revel et al,. (2019,2021)]
# ==============================================================================================
# References:
# 1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating 
# global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 
# 1–34. https://doi.org/10.1029/2020wr027876
# 2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical 
# Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: 
# A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829
# 3. Revel, M., Zhou, X., Yamazaki, D., & Kanae, S. (2023). Assimilation of transformed water 
# surface elevation to improve river discharge estimation in a continental-scale river. 
# Hydrology and Earth System Sciences, 27(3), 647–671. https://doi.org/10.5194/hess-27-647-2023
# ==============================================================================================
# created by Ikeshima & Menaka
# Menaka@IIS 2021
#***********************************************************************************************

################################################################################################
#
# This program runs the whole data assimilation in mool server@IIS.
#
# Please change the necessary PBS requirments below
#
################################################################################################

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F40S
#PBS -l select=1:ncpus=40:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N HydroDA

#source ~/.bashrc
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

# get number of cpus
#export NCPUS=`cat ${PBS_NODEFILE} | wc -l`
NCPUS=40

# OMP Settings
export OMP_NUM_THREADS=$NCPUS

# go to working dirctory
HydroDA="/cluster/data6/menaka/HydroDA"
# HydroDAout="/cluster/data7/menaka/HydroDA"
HydroDAout="/cluster/data8/menaka/HydroDA"
mkdir -p $HydroDAout
# HydroDAout="/work/a06/menaka/HydroDA"

#----------
# cd $HydroDA
cd $HydroDAout
#cd $PBS_O_WORKDIR
#cd $swotda

#******************************************************************************************
# experiment : edit the experiment name in here it will be written in $HydroDA/$EXP/exp.txt
# before running run_mool.sh , please edit the necessary experimental settings in params.py
#******************************************************************************************

#====================================================================
# experiment name [XXX_YYY_ZZZ_WWW_NNN]
# 1. Assimilation method [direct(DIR), anomaly(ANO), normalized(NOM)]
# 2. Observation variable [WSE/DIS]
# 3. Runoff Data [e.g., E2O, VICBC, ECMWF]
# 4. Observation data [e.g., HydroWeb(HWEB), CGLS] 
# 5. Number for identifying the experiment [e.g., 001]: 0XX - regional, 1XX - global
#====================================================================
for expnum in `seq 19 20`; # 2023/7/10 @9:00 E20
do
    num1=$(($expnum + 50))
    num1=$(printf "%03g" $num1)

    EXP="ANO_WSE_ISIMIP3a_SWOT_$num1" # for SWOTH08 
    # name refernce
    # 1 - no parameter error
    # 2 - rivght error
    # 3 - rivwth error
    # 4 - rivman error
    # 5 - fldhgt error
    # 6 - all paremeter error
    # from 51-70 - 20 multiple true experiments

    echo $EXP, `pwd`

    cd $HydroDAout"/out/"$EXP

    # run the main code using virtual environment
    # run main code
    touch ./__init__.py
    python run.py
done

wait

conda deactivate