#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import sys
from sys import *

argvs = sys.argv
year=int(argvs[1])
month=int(argvs[2])
date=int(argvs[3])
obj=argvs[4]

#calc leap year
leap_y=0 #0:not 1:leap year
if year%4==0:
    if year%100==0 and year%400!=0:
        leap_y=0
    else:
        leap_y=1

#calc last day of the month
if month==1:
    last_d=31
if month==2:
    if leap_y==1:
        last_d=29
    else:
        last_d=28
if month==3:
    last_d=31
if month==4:
    last_d=30
if month==5:
    last_d=31
if month==6:
    last_d=30
if month==7:
    last_d=31
if month==8:
    last_d=31
if month==9:
    last_d=30
if month==10:
    last_d=31
if month==11:
    last_d=30
if month==12:
    last_d=31

next_y=year
next_m=month
next_d=date
#calc next date
if date+1>last_d:
    next_d=1
    if month+1>12:
        next_m=1
        next_y=year+1
    else:
        next_m=month+1
else:
    next_d=date+1

if obj=="year":
    print next_y
if obj=="month":
    print next_m
if  obj=="date":
    print next_d
