#/bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=60gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N HydroDA

#source ~/.bashrc

#export OMP_NUM_THREADS=40

# got to working dirctory
HydroDA="/cluster/data6/menaka/HydroDA"

#--
cd $HydroDA
#cd $PBS_O_WORKDIR
#cd $swotda

# experiment : edit the experiment name in here and params.py experiment()
# before running run_mool.sh , please edit the nessary experimental settings in params.py
EXP="E2O_HydroWeb16"
#IFACTOR="1.08"

mkdir $HydroDA"/out"
mkdir $HydroDA"/out/"$EXP

#write experiment name
echo $EXP > $HydroDA"/out/"$EXP"/exp.txt"

# copy params.py
cp -r "$HydroDA/gosh/params.py" "$HydroDA/out/$EXP/params.py"
#cp -r "$HydroDA/gosh/params.py" "$HydroDA/src/params.py"

cp -r $HydroDA"/src/"run.py $HydroDA"/out/"$EXP"/"run.py
ln -sf $HydroDA"/src/"main_code.py $HydroDA"/out/"$EXP"/"main_code.py

cd $HydroDA"/out/"$EXP
#python run.py
# run the main code using virtual environment
/home/menaka/miniconda3/envs/pydef/bin/python2.7 run.py
