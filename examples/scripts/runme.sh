#! /bin/bash 
export OMP_NUM_THREADS=4 
for Temp in 010 100 200 300 400 500 600 700 800 900 950 1000 1050 1100 1150 1200 1250 1300 1500
do
 mkdir T$Temp/ 2>/dev/null
 echo "Temp: " $Temp
 cp Base/* T$Temp/
 cd T$Temp/
 sed -i "s/TEMP/$Temp/g" inpsd.dat
 if [ $Temp -ne 010 ]
 then
   sed -i "s/#restart/restart/" inpsd.dat 
   cp ../restart.bccFe100.out .
 fi
 ../../../source/sd > out.log
 cp restart.bccFe100.out ../
 cd ..
done
exit
