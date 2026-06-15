#!/bin/sh

clean=$1  #  argument for cleaning [yes/clean]

#HydroDA="cluster/data6/menaka/HydroDA"
HydroDA=`pwd`/..
BASE=$HydroDA  # code location

libs="img src"
for lib in $libs
do
  cd $BASE/$lib
  if [ $clean = "yes" ];then
    make clean
    echo "*********** $lib **********"
    make -B all
  elif [ $clean = "clean" ];then
    make clean
  else
    echo "*********** $lib **********"
    make all
  fi
done
