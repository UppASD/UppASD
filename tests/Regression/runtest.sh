#!/bin/bash
rm -f *out
export OMP_NUM_THREADS=1
$SD_BINARY | tee out
unset OMP_NUM_THREADS
