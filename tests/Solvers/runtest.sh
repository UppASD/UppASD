#!/bin/bash

for solver in 1 4 5
do
  sed "s/SOLVER/$solver/g" inpsd.dat.base > inpsd.dat
  echo $solver
  $SD_BINARY > out.$solver
  #/home/andersb/SD/UppASD_3.2/source/sd > out.$solver
  cp averages.SolvTest.out averages.$solver 2> /dev/null
  cp totenergy.SolvTest.out totenergy.$solver 2> /dev/null
rm *.out meminfo tempfile 2> /dev/null
done
rm -f inpsd.dat
