#! /bin/bash 
for Temp in 100 200 300 400 500 600 700 800 900 950 1000 1100 1200 1300 1400 1500 1600
do
    mkdir T$Temp/ 2>/dev/null
    echo "Temp: " $Temp
    cp Base/* T$Temp/
    cd T$Temp/
    sed -i "s/TEMP/$Temp/g" inpsd.dat
    if [ $Temp -ne 100 ]
    then
	sed -i "s/#restart/restart/" inpsd.dat 
	cp ../restart.bccFe100.out .
    fi
    ${SD_BINARY} > out.log
    cp restart.bccFe100.out ../
    cd ..
done
exit
