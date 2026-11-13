#!/bin/bash
export OMP_NUM_THREADS=1
$SD_BINARY | tee out
