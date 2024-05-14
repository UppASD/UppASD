#!/bin/bash

for n in 001 002 003

do

    mkdir disp-$n
    cd disp-$n
    cp -p ../POSCAR-$n POSCAR
    cp -p ../INCAR .
    cp -p ../KPOINTS .
    cp -p ../POTCAR .
    cp -p ../jobvaspN4n128.sh .
    sbatch jobvaspN4n128.sh
    cd ..
    
done
