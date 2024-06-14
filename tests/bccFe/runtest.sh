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
   if [ "$mode" != "S" ]
   then
       sed -i'' -e "s/ip_nphase/#ip_nphase/g" inpsd.dat
       sed -i'' -e "s/Nstep/mcNstep/g" inpsd.dat
   else
       sed -i'' -e "s/ip_mcanneal/#ip_mcanneal/g" inpsd.dat
   fi
   $SD_BINARY | tee out.$mode
   cp ./averages.bcc_Fe_T.out averages.$mode
   cp ./cumulants.bcc_Fe_T.out cumulants.$mode
   rm *.out *.json *.yaml
done

rm -f inpsd.dat
