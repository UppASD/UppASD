#!/bin/bash

for mode in S M H
do
   #if [ "$mode" != "S" ]
   #then
   #    export OMP_NUM_THREADS=1
   #else
   #    unset OMP_NUM_THREADS
   #fi
   #export OMP_NUM_THREADS=1

   export OMP_NUM_THREADS=1
   sed "s/MODE/$mode/g" inpsd.dat.base > inpsd.dat
   #/home/andersb/SD/UppASD_3.2/source/sd | tee out.$mode
   $SD_BINARY | tee out.$mode
   cp ./averages.bcc_Fe_T.out averages.$mode
   cp ./cumulants.bcc_Fe_T.out cumulants.$mode
   rm *.out
done

rm -f inpsd.dat
