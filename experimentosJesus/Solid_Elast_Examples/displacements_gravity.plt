reset
#
#set term postscript enhanced color
#set output "deslocamentos_problm_grav.eps"
set xlabel 'Y Axis'
set ylabel 'U_y'
set grid
set key left
set title "Vertical Displacements Solutions for Body subject to Gravity"
plot "dxelast01.ht01/displac.dat" using ($0/2):($2) every 26 title "Numerical" with linespoints, "dxelast01.ht01/fnodes.stoc" using ($2):(($2)*($2)-200*200)/(2000.0) every 26 title "Analytical" with lines lw 2 
