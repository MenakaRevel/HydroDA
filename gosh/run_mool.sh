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

source ~/.bashrc

export OMP_NUM_THREADS=20

# got to working dirctory
SWOTDA="/cluster/data6/menaka/SWOTDA_sensivity"

#--
cd $SWOTDA
#cd $PBS_O_WORKDIR
#cd $swotda

# experiment
EXP="ECMWF_womc_biased_if_fixed1.08"
IFACTOR="1.08"

# run the main code using virtual environment
/home/menaka/miniconda3/envs/pydef/bin/python2.7 PBS_run.py $EXP $IFACTOR
