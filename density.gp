
plot 'solution' u 1:2 t 'density' w p, 'riemann.data' u 1:2 t 'exact' w l
set term postscript eps enhanced mono
set output 'density.eps'
replot
set output
set terminal x11
