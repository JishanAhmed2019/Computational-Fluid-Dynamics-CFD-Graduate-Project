#set xrang [0.0:1.0]
#set yrang [0.0:1.0]
plot 'solution' u 1:4 t 'pressure' w p, 'riemann.data' u 1:4 t 'exact' w l
set term postscript eps enhanced mono
set output 'pressure.eps'
replot
set output
set terminal x11
