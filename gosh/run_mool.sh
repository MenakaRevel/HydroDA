#!/bin/bash
#***********************************************************************************************
# Data Assimilation using LETKF and empircal local patches [Revel et al,. (2019,2021)]
# ==============================================================================================
# Reference:
# 1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating 
# global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 
# 1–34. https://doi.org/10.1029/2020wr027876
# 2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical 
# Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: 
# A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829
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
#PBS -q E40
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
NCPUS=48

# OMP Settings
export OMP_NUM_THREADS=$NCPUS

# go to working dirctory
HydroDA="/cluster/data6/menaka/HydroDA"
HydroDAout="/cluster/data7/menaka/HydroDA"

#----------
# cd $HydroDA
cd $HydroDAout
#cd $PBS_O_WORKDIR
#cd $swotda

#******************************************************************************************
# experiment : edit the experiment name in here it will be written in $HydroDA/$EXP/exp.txt
# before running run_mool.sh , please edit the necessary experimental settings in params.py
#******************************************************************************************

# experiment name
# EXP="VIC_BC_HydroWeb11"
# EXP="E2O_HydroWeb23"
# EXP="test_wse"

#====================================================================
# experiment name [XXX_YYY_ZZZ_WWW]
# 1. Assimilation method [direct(DIR), anomaly(ANO), normalized(NOM)]
# 2. Observation variable [WSE/DIS]
# 3. Runoff Data [e.g., E2O, VICBC, ECMWF]
# 4. Observation data [HydroWeb(HWEB)] 
# 5. Number for identifying the experiment [e.g., 001]: 0XX - regional, 1XX - global

# EXP="DIR_WSE_ECMWF_HWEB_014"
# EXP="ANO_WSE_ECMWF_HWEB_012"
# EXP="NOM_WSE_ECMWF_HWEB_012"
# EXP="test_virtual"
# EXP="test_wse"
EXP="NOM_WSE_E2O_HWEB_101" # for glb_15min
EXP="DIR_WSE_E2O_SWOT_001" # for SWOTH08

# mkdir -p $HydroDA"/out/"$EXP
mkdir -p $HydroDAout"/out/"$EXP

# go to working directory
# cd $HydroDA"/out/"$EXP
cd $HydroDAout"/out/"$EXP

#write experiment name
# echo $EXP > $HydroDA"/out/"$EXP"/exp.txt"
echo $EXP > "./exp.txt"

#write NCPUS
# echo $NCPUS > $HydroDA"/out/"$EXP"/ncpus.txt"
echo $NCPUS > $"./ncpus.txt"

# copy params.py
# cp -r $HydroDA/gosh/params.py     $HydroDA/out/$EXP/params.py
# cp -r $HydroDA/gosh/params_virt.py     $HydroDA/out/$EXP/params.py
# cp -r $HydroDA/gosh/params_real.py     $HydroDA/out/$EXP/params.py
# cp -r $HydroDA/gosh/params_real.py     ./params.py
cp -r $HydroDA/gosh/params_virt.py     ./params.py

# copy running related files
# cp -r $HydroDA/src/run.py           $HydroDA/out/$EXP/run.py
# cp -r $HydroDA/src/main_code.py     $HydroDA/out/$EXP/main_code.py
# cp -r $HydroDA/src/prep_init.py     $HydroDA/out/$EXP/prep_init.py
# cp -r $HydroDA/src/prep_runoff.py   $HydroDA/out/$EXP/prep_runoff.py
# cp -r $HydroDA/src/prep_obs.py      $HydroDA/out/$EXP/prep_obs.py
# cp -r $HydroDA/src/wrt_expset.py    $HydroDA/out/$EXP/wrt_expset.py

cp -r $HydroDA/src/run.py           ./run.py
cp -r $HydroDA/src/main_code.py     ./main_code.py
cp -r $HydroDA/src/prep_init.py     ./prep_init.py
cp -r $HydroDA/src/prep_runoff.py   ./prep_runoff.py
cp -r $HydroDA/src/prep_obs.py      ./prep_obs.py
cp -r $HydroDA/src/wrt_expset.py    ./wrt_expset.py


# run the main code using virtual environment
# run main code
touch ./__init__.py &
python run.py &

wait

conda deactivate