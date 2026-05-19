#!/bin/bash

for Temp in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0

do
    mkdir T$Temp/
    echo "Temp: " $Temp
    cp Base/* T$Temp/
    cd T$Temp/
    sed -i "s/TEMP/$Temp/g" inpsd.dat
    ${SD_BINARY} > out.log
    cd ..
done
