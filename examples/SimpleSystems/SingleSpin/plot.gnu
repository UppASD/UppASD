#! /usr/bin/gnuplot

set terminal pngcairo enhanced size 2560,1600
set output 'GilbertTotal.png'
unset grid
unset key
set style line 1 lt -1 lc rgb "#7fc97f" lw 2.5 pt 1 ps 1.25 dashtype 0
set style line 2 lt  0 lc rgb "#beaed4" lw 2.5 pt 4 ps 1.25 dashtype 1 
set style line 3 lt  1 lc rgb "#fdc086" lw 2.5 pt 6 ps 1.25 dashtype 2 
set style line 4 lt  2 lc rgb "#386cb0" lw 2.5 pt 8 ps 1.25 dashtype 3 
set style line 5 lt  3 lc rgb "#f0027f" lw 2.5 pt 10 ps 1.25 dashtype 4 
set style line 6 lt  4 lc rgb "#bf5b17" lw 2.5 pt 12 ps 1.25 dashtype 5 
set style line 7 lt  5 lc rgb "#666666" lw 2.5 pt 27 ps 1.25 dashtype 0 
set multiplot title "\n Magnetization Vs. Time for different dampings" font ", 40"
set tics font ",24"

set border 15 lw 5
set xlabel "Time [ns]" font ", 24"
unset title
set size 0.5, 0.45
set origin 0.0, 0.45
set ylabel "M_x" font ",24"
plot "./DAMP0.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 1 ti "{/Symbol a}=0","./DAMP0.0100/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 3 ti "{/Symbol a}=0.01","./DAMP1.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 5 ti "{/Symbol a}=1","./DAMP100.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 7 ti "{/Symbol a}=100" #"./DAMP0.0010/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 2 ti "{/Symbol a}=0.001", "./DAMP0.0100/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 3 ti "{/Symbol a}=0.01", "./DAMP0.1000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 4 ti "{/Symbol a}=0.1", "./DAMP1.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 5 ti "{/Symbol a}=1", "./DAMP10.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 6 ti "{/Symbol a}=10","./DAMP100.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):2 w lp ls 7 ti "{/Symbol a}=100"
unset ylabel

unset title
set size 0.5, 0.45
set origin 0.5, 0.45
set ylabel "M_y" font ", 24"
plot "./DAMP0.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 1 ti "{/Symbol a}=0","./DAMP0.0100/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 3 ti "{/Symbol a}=0.01","./DAMP1.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 5 ti "{/Symbol a}=1","./DAMP100.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 7 ti "{/Symbol a}=100" #"./DAMP0.0010/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 2 ti "{/Symbol a}=0.001", "./DAMP0.0100/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 3 ti "{/Symbol a}=0.01", "./DAMP0.1000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 4 ti "{/Symbol a}=0.1", "./DAMP1.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 5 ti "{/Symbol a}=1", "./DAMP10.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 6 ti "{/Symbol a}=10","./DAMP100.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):3 w lp ls 7 ti "{/Symbol a}=100"
unset ylabel


unset title
set size 0.6, 0.45
set origin 0.3, 0.005
set ylabel "M_z" font ", 24"
set key outside font ",24"
set yrange [0.49:1.01]
plot "./DAMP0.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 1 ti "{/Symbol a}=0","./DAMP0.0100/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 3 ti "{/Symbol a}=0.01","./DAMP1.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 5 ti "{/Symbol a}=1","./DAMP100.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 7 ti "{/Symbol a}=100" #"./DAMP0.0010/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 2 ti "{/Symbol a}=0.001", "./DAMP0.0100/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 3 ti "{/Symbol a}=0.01", "./DAMP0.1000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 4 ti "{/Symbol a}=0.1", "./DAMP1.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 5 ti "{/Symbol a}=1", "./DAMP10.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 6 ti "{/Symbol a}=10","./DAMP100.0000/averages.SingleSP.out" u ($1*1e-14/1e-9):4 w lp ls 7 ti "{/Symbol a}=100"
unset ylabel
