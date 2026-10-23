#! /bin/bash
m_max=`tail -n1 MT*/cumulants.*.out | grep E | awk '{ print $2 }' | sort -gr | head -1 `
cv_max=`tail -n1 MT*/cumulants.*.out | grep E | awk '{ print $7 }' | sort -gr | head -1 `
#cv_max=1
x_max=`tail -n1 MT*/cumulants.*.out | grep E | awk '{ print $6 }' | sort -gr | head -1 `
echo "# Temp.   Mavg     UBinder    Susc.      Cv" > thermal.dat
echo "# Temp.   Mavg     UBinder    Susc.      Cv" > thermal.norm.dat
for Temp in 001 010 100 200 300 400 500 600 700 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 1200 1300 1400 1500 1600
 do
          tail -1 MT$Temp/cumulants.*.out | awk '{ printf "%6i  %f  %f  %f  %f\n", temp, $2,$5,$6,$7 }' temp=$Temp >> thermal.dat
          tail -1 MT$Temp/cumulants.*.out | awk '{ printf "%6i  %f  %f  %f  %f\n", temp, $2/m_max,$5,$6/x_max,$7/cv_max }' m_max=$m_max x_max=$x_max cv_max=$cv_max temp=$Temp >> thermal.norm.dat
 done
exit
