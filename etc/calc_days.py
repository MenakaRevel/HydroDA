#!/opt/local/bin/python
# -*- coding: utf-8 -*-
# import numpy as np
import datetime
import sys
import os
#--
syear=int(sys.argv[1])
smon =int(sys.argv[2])
sday =int(sys.argv[3])
eyear=int(sys.argv[4])
emon =int(sys.argv[5])
eday =int(sys.argv[6])
#--
start_dt=datetime.date(syear,smon,sday)
last_dt=datetime.date(eyear,emon,eday)
last=(last_dt-start_dt).days + 1
#**********
N = int(last)
print (N)