#! /usr/bin/gnuplot

set terminal pngcairo enhanced size 2560,1600
set output 'relaxation.png'
unset grid
set style line 1 lt -1 lc rgb "#7fc97f" lw 5 pt 1 ps 2 dashtype 0
set style line 2 lt  0 lc rgb "#beaed4" lw 5 pt 4 ps 2 dashtype 1 
set style line 3 lt  1 lc rgb "#fdc086" lw 5 pt 6 ps 2 dashtype 2 
set style line 4 lt  2 lc rgb "#386cb0" lw 5 pt 8 ps 2 dashtype 3 
set style line 5 lt  3 lc rgb "#f0027f" lw 5 pt 10 ps 2 dashtype 4 
set style line 6 lt  4 lc rgb "#bf5b17" lw 5 pt 12 ps 2 dashtype 5 
set style line 7 lt  5 lc rgb "#666666" lw 5 pt 27 ps 2 dashtype 0 
set tics font ",30"

set size 0.9,0.9
set origin 0.05, 0.05
set border 15 lw 5
set xlabel "Time [ns]" font ", 30"
set logscale x
unset title
set ylabel "M_z" font ",30" offset "-5,0"
set key top right font ",30"
plot "./T1/EN200/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 1 ti "T=1 K", "./T2/EN200/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 2 ti "T=2 K", "./T3/EN200/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 3 ti "T=3 K", "./T4/EN200/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 4 ti "T= 4 K", "./T5/EN200/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 5 ti "T=5 K","./T10/EN200/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 6 ti "T=10 K" 

reset
set terminal pngcairo enhanced size 2560,1600

set output 'singleTemp.png'
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
set ylabel "M_z" font ",30" offset "-5,0"

p "./T4/EN1/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 4 noti
