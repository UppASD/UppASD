unset key
unset colorbox
#unset xtics
#unset ytics
unset border
set xlabel "qx"
set ylabel "qy"
#set palette rgbformula  6,6,6
set palette rgbformula  -6,-6,-6
#set palette rgbformula  30,31,32
#set palette rgbformula  33,13,10
set pm3d map
#set zrange [0.00:1]
set zrange [-0:1]
set cbrange [-0:1]
#splot "sqw.mat" matrix w pm3d
splot "sqw.mat" matrix using ($1*2*pi/3.16/20):($2*YSCALE):($3) w pm3d
pause -1
set term png transparent
set output "sqw.png"
splot "sqw.mat" matrix using ($1*2*pi/3.16/20):($2*YSCALE):($3) w pm3d

