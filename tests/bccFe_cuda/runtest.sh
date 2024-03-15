#!/bin/bash
SD_BINARY=../../bin/sd.f95.cuda
export OMP_NUM_THREADS=1
$SD_BINARY | tee out.test
cp ./averages.bcc_Fe_T.out averages.test
cp ./cumulants.bcc_Fe_T.out cumulants.test
rm -f *.out *.json *.yaml

