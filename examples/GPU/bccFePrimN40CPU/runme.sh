#! /bin/bash
for Temp in 001 010 100 200 300 400 500 600 700 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 1200 1300 1400 1500 1600
do
    # Initial phase
    mkdir IT$Temp/ 2>/dev/null
    echo "Temp: " $Temp
    cp Base/* IT$Temp/
    cd IT$Temp/
    sed -i "s/NX/40/g" inpsd.dat
    sed -i "s/NY/40/g" inpsd.dat
    sed -i "s/NZ/40/g" inpsd.dat
    sed -i "s/NSTEP/10000/g" inpsd.dat
    sed -i "s/TEMP/$Temp/g" inpsd.dat
    sed -i "s/DAMP/1.00/g" inpsd.dat
    if [ $Temp -ne 001 ]
    then
	sed -i "s/INITMAG/4/" inpsd.dat
	sed -i "s/#restart/restart/" inpsd.dat
	cp ../restart.bccFe000.out .
    else
	sed -i "s/INITMAG/3/" inpsd.dat
    fi
    ${SD_BINARY} > out.log
    cp restart.bccFe000.out ../
    cd ..

    # Measurement phase
    mkdir MT$Temp/ 2>/dev/null
    echo "Temp: " $Temp
    cp Base/* MT$Temp/
    cd MT$Temp/
    sed -i "s/NX/40/g" inpsd.dat
    sed -i "s/NY/40/g" inpsd.dat
    sed -i "s/NZ/40/g" inpsd.dat
    sed -i "s/NSTEP/100000/g" inpsd.dat
    sed -i "s/TEMP/$Temp/g" inpsd.dat
    sed -i "s/DAMP/1.00/g" inpsd.dat
    sed -i "s/INITMAG/4/" inpsd.dat
    sed -i "s/#restart/restart/" inpsd.dat
    cp ../restart.bccFe000.out .
    ${SD_BINARY} > out.log
    cp restart.bccFe000.out ../
    cd ..
done
exit
