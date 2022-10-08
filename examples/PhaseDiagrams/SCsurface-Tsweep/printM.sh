#! /bin/csh -f 
foreach Temp (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)
          echo $Temp | sed ' s/K//' > tmp
          cat tmp >> Tlist
          cat T$Temp/cumulants.scHeis64.out | tail -1 | awk '{print $2}' > tmp1
          cat T$Temp/cumulants.scHeis64.out | tail -1 | awk '{print $5}' > tmp2
          cat tmp1 >> Mlist
          cat tmp2 >> Ulist
end

paste Tlist Mlist Ulist > magnetization.dat
rm Tlist Mlist Ulist tmp tmp1 tmp2
exit
