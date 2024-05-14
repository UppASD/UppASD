#!/bin/bash

for n in 001

do

    mkdir disp-$n
    cd disp-$n
    cp -p ../POSCAR-$n POSCAR
    cp -p ../INCAR .
    cp -p ../KPOINTS .
    cp -p ../POTCAR .
    cp -p ../jobvaspN4n16t8.sh .
    sbatch jobvaspN4n16t8.sh
    cd ..
    
done
