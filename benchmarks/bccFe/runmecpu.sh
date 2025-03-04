#! /bin/bash
for nx in 10 20 30 40 50 60
do

    # ASD CPU
    mkdir ASDCPUN$nx/ 2>/dev/null
    echo "NX: " $nx
    cp Base/* ASDCPUN$nx/
    cd ASDCPUN$nx
    sed -i "s/NX/$nx/g" inpsd.dat
    sed -i "s/NY/$nx/g" inpsd.dat
    sed -i "s/NZ/$nx/g" inpsd.dat
    sed -i "s/MODE/S/g" inpsd.dat
    sed -i "s/GPU/0/g" inpsd.dat
    ../../../bin/sd.f95.cuda > out.log
    cd ..

    # MC CPU
    mkdir MCCPUN$nx/ 2>/dev/null
    echo "NX: " $nx
    cp Base/* MCCPUN$nx/
    cd MCCPUN$nx
    sed -i "s/NX/$nx/g" inpsd.dat
    sed -i "s/NY/$nx/g" inpsd.dat
    sed -i "s/NZ/$nx/g" inpsd.dat
    sed -i "s/MODE/S/g" inpsd.dat
    sed -i "s/nsteps/mcnsteps/g" inpsd.dat
    sed -i "s/GPU/0/g" inpsd.dat
    ../../../bin/sd.f95.cuda > out.log
    cd ..

done
exit
