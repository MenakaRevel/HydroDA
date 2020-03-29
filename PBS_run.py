#!/opt/local/bin/python
# -*- coding: utf-8 -*-

########################
#
# this program run the whole program
#
########################

# run letkf.py
# 2018-12-29 @ Menaka
import sys
import main_code

exp=sys.argv[1]
ifact=sys.argv[2]

try:
  main_code.main_act(exp,ifact)
except Exception as e:
  print e


## check before run ################
#
# 1. compile CaMa-Flood
# 2. set CaMa/src/control_inp.F propertly
# 3. set parameters
# 4. set local patches and the weights
# 4. compile data_assim.f90 using by executing compileMKL.sh at from this directory
#
#
