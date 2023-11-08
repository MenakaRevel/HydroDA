#!/bin/bash
### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N prep_obs
#==
USER=`whoami`

# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

# OMP Settings
NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/HydroDA/etc"

python prep_SWOT_obs.py