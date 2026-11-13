#!/bin/bash

export OMP_NUM_THREADS=1
$SD_BINARY | tee out.log
cp ./averages.FeCo__B2.out averages.test
cp ./cumulants.FeCo__B2.out cumulants.test
rm -f *.out *.json *.yaml

