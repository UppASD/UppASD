#!/bin/bash
#/home/andersb/SD/UppASD_3.2/source/sd | tee out
export OMP_NUM_THREADS=1
$SD_BINARY | tee out
