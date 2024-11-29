#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#*************************************************************************************
# Data Assimilation using LETKF and empircal local patches [Revel et al,. (2019,2021)]
# ====================================================================================
# Reference:
# 1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating 
# global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 
# 1–34. https://doi.org/10.1029/2020wr027876
# 2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical 
# Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: 
# A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829
# ====================================================================================
# created by Ikeshima & Menaka
# Menaka@IIS 2021
#*************************************************************************************

########################################
#
# This program runs the whole algorithm
#
########################################

## check before run ################
#
# 1. compile CaMa-Flood
# 2. set parameters params.py@gosh directory
# 3. set local patches folder
# 4. compile fortran codes in ./src folder %sh compile.sh "yes"

import sys
import main_code
import prep_init as init 
import prep_runoff as inpt
import prep_obs as obs
import wrt_expset as expset

try:
  # write the experimental settings in log file
  print ("write experiment log file")
  expset.write_text()

  # make necessary directories
  print ("initializing....")
  init.initial()

  # prepare runoff ensembles
  print ("prepare input")
  inpt.prepare_input()

  # prepare observations
  print ("prepare observation")
  obs.prepare_obs()

  # initial inflation parameter rho for assimilation
  print ("make intial inflation")
  init.make_initial_infl()

  # prepare the mean and std for anomaly/normalized assimilation
  print ("save statistics")
  init.save_statistic()

  # run main code
  print ("running main assimilation code....")
  main_code.main_act()

  print ("HydroDA completed.")
except Exception as e:
  print (e)