
plot 'solution' u 1:3  t 'velocity' w p,  'riemann.data' u 1:3 t 'exact' w l
set term postscript eps enhanced mono
set output 'velocity.eps'
replot
set output
set terminal x11
