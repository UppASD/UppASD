#! /bin/csh -f 

foreach Temp (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)

mkdir T$Temp/
echo "Temp: " $Temp
#T=`echo $Temp | awk '{ printf "%.4i\n",$1 }' `
cp Base/* T$Temp/
cd T$Temp/
sed -i "s/TEMP/$Temp/g" inpsd.dat
#sbatch runme.sh
../../../../source/sd
cd ..

end
exit
