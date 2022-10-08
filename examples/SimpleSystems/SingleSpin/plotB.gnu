#! /usr/bin/gnuplot

set terminal pngcairo enhanced size 2560,1600
set output 'Switching.png'
unset grid
set style line 1 lt -1 lc rgb "#7fc97f" lw 3 pt 1 ps 2 dashtype 0
set style line 2 lt  0 lc rgb "#beaed4" lw 3 pt 4 ps 2 dashtype 1 
set style line 3 lt  1 lc rgb "#fdc086" lw 3 pt 6 ps 2 dashtype 2 
set style line 4 lt  2 lc rgb "#386cb0" lw 3 pt 8 ps 2 dashtype 3 
set style line 5 lt  3 lc rgb "#f0027f" lw 3 pt 10 ps 2 dashtype 4 
set style line 6 lt  4 lc rgb "#bf5b17" lw 3 pt 12 ps 2 dashtype 5 
set style line 7 lt  5 lc rgb "#666666" lw 3 pt 27 ps 2 dashtype 0 
set tics font ",30"

set size 0.9,0.9
set origin 0.05, 0.05
set border 15 lw 5
set xlabel "Time [ns]" font ", 30"
unset title
set yrange [-1.05:1.05]
set ylabel "M_z" font ",30"
set key center right font ",30"
plot "./FIELD0.00/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 1 ti "B=0.00 T", "./FIELD0.20/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 2 ti "B=0.20 T", "./FIELD0.25/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 3 ti "B=0.25 T", "./FIELD0.35/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 4 ti "B=0.35 T", "./FIELD0.45/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 5 ti "B=0.45 T" #,"./FIELD0.30/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 4 ti "B=0.30 T", "./FIELD0.40/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 4 ti "B=0.40 T"

