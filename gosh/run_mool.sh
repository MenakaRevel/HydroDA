#!/bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=40gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N SWOTDA_1.08

#source ~/.bashrc

export OMP_NUM_THREADS=20

# got to working dirctory
HydroDA="/cluster/data6/menaka/HydroDA"

#--
cd $HydroDA
#cd $PBS_O_WORKDIR
#cd $swotda

# experiment
EXP="E2O_womc"
IFACTOR="1.08"

mkdir $HydroDA"/out"
mkdir $HydroDA"/out/"$EXP

# copy params.py 
cp -r "$HydroDA/gosh/params.py" "$HydroDA/out/$EXP/params.py"
cp -r "$HydroDA/gosh/params.py" "$HydroDA/src/params.py"

ln -sf $HydroDA"/src/"run.py $HydroDA"/out/"$EXP"/"run.py
#ln -sf $HydroDA"/src/"main_code.py $HydroDA"/out/"$EXP"/"main_code.py
# run the main code using virtual environment
#/home/menaka/miniconda3/envs/pydef/bin/python2.7 $HydroDA"/src/"PBS_run.py $EXP $IFACTOR
cd $HydroDA"/out/"$EXP
python run.py
