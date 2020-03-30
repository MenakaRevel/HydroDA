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
# 2. set parameters
# 3. set local patches folder
# 4. compile fortran codes in ./src folder
#
#
