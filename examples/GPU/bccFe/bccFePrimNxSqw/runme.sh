#!/bin/bash
for Nx in 40 48 60 100
do
for Temp in 001 #010 100 300 1000
do
    # Sampling of S(q,w) at Gilbert damping a=0.001
    mkdir N$Nx''MT$Temp 2>/dev/null
    echo "Cell size: " $Nx "Temp: " $Temp
    cp Base/* N$Nx''MT$Temp
    cd N$Nx''MT$Temp
    sed -i "s/NX/$Nx/g" inpsd.dat
    sed -i "s/NY/$Nx/g" inpsd.dat
    sed -i "s/NZ/$Nx/g" inpsd.dat
    sed -i "s/TEMP/$Temp/g" inpsd.dat
    sed -i "s/IPNSTEP/10000/g" inpsd.dat
    sed -i "s/IPDAMP/1.0/g" inpsd.dat
    sed -i "s/NSTEP/20000/g" inpsd.dat
    sed -i "s/DAMP/0.001/g" inpsd.dat
    ${SD_BINARY} > out.log
    cd ..
done
done
exit
