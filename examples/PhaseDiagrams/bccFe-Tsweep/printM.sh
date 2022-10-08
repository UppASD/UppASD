#! /bin/bash
m_max=`tail -n1 T*/cumulants.*.out | grep E | awk '{ print $2 }' | sort -gr | head -1 `
cv_max=`tail -n1 T*/cumulants.*.out | grep E | awk '{ print $7 }' | sort -gr | head -1 `
#cv_max=1
x_max=`tail -n1 T*/cumulants.*.out | grep E | awk '{ print $6 }' | sort -gr | head -1 `
echo "# Temp.   Mavg     UBinder    Susc.      Cv" > thermal.dat
echo "# Temp.   Mavg     UBinder    Susc.      Cv" > thermal.norm.dat
for Temp in 100 200 300 400 500 600 700 800 900 950 1000 1100 1200 1300 1400 1500 1600
 do
          tail -1 T$Temp/cumulants.*.out | awk '{ printf "%6i  %f  %f  %f  %f\n", temp, $2,$5,$6,$7 }' temp=$Temp >> thermal.dat
          tail -1 T$Temp/cumulants.*.out | awk '{ printf "%6i  %f  %f  %f  %f\n", temp, $2/m_max,$5,$6/x_max,$7/cv_max }' m_max=$m_max x_max=$x_max cv_max=$cv_max temp=$Temp >> thermal.norm.dat
 done
exit
