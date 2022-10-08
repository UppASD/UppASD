set xlabel "Temperature (K)"
set xrange [0:1600]
set yrange [0:1.1]
p "thermal.norm.dat" u 1:2 w lp ti "Magnetization"
rep "thermal.norm.dat" u 1:3 w lp ti "Binder cumulant"
rep "thermal.norm.dat" u 1:4 w lp ti "Susceptibility"
rep "thermal.norm.dat" u 1:5 w lp ti "Specific heat"
pause(-1)
set term png
set output "bccfe.png"
rep
