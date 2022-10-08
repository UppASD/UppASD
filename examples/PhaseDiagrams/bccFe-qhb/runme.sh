#! /bin/bash 
echo "Calculating adiabatic magnon DOS"
mkdir AMS
cp Base/* AMS
cd AMS
mv inpsd.dat.ams inpsd.dat
../../../../source/sd > out.log
cd ..

echo "Performing temperature sweep"
for Temp in 001 100 200 300 400 500 600 700 800 900 950 1000 1100 1200 1300 1400 1500
do
    mkdir T$Temp/ 2>/dev/null
    echo "Temp: " $Temp
    cp Base/* T$Temp/
    cd T$Temp/
    cp ../AMS/magdos.bccFe100.out magdos.bccFe.dat
    sed "s/TEMP/$Temp/g" inpsd.dat.mc > inpsd.dat
    ../../../../source/sd > out.log
    cp restart.bccFe100.out ../
    cd ..
done
exit
