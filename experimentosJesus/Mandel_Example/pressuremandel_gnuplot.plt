reset
#
## Define configurations for graphics display
#
set grid
set zeroaxis lt -1 lw 3
#set pointsize 2


set yrange [-0.1:0.6]
set xrange [-0.1:1.1]
set ylabel "Pressure"
set xlabel "X Axis"


set key outside 

#
## Define output parameters
#set terminal postscript eps enhanced color solid
#set terminal png 
#set terminal postscript monochrome dashed
#
#set output 'step1.ps'
#
F = 1.0e+8
A = 1.0e+2
B = 1.0e+1
xnorm = A/F
#
## Plot de functions
#
## Set Title
#
set title "Pressure Solutions for Mandel Problem \n Time 1.0E3"
#set output 'mandeltime_1x10E3.eps'
plot "dxelast10.ht01/prsr000.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Init Pressr' with lines lt -1 , "dxelast10.ht01/prsr001.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Numerical' with linespoints lt 1, "analytical/mandelteor10.3dat" using ($1):($2) with lines lt 2 title 'Analytical'
pause -1
#set output 'mandeltime_1x10E4.eps'
set title "Pressure Solutions for Mandel Problem \n Time 1.0E4"
plot "dxelast10.ht01/prsr000.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Init Pressr' with lines lt -1 ,"dxelast10.ht01/prsr010.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Numerical' with linespoints lt 1, "analytical/mandelteor10.4dat" using ($1):($2) with lines lt 2 title 'Analytical'
pause -1
#set output 'mandeltime_5x10E4.eps'
set title "Pressure Solutions for Mandel Problem \n Time 5.0E4"
plot "dxelast10.ht01/prsr000.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Init Pressr' with lines lt -1 ,"dxelast10.ht01/prsr050.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Numerical' with linespoints lt 1, "analytical/mandelteor510.4dat" using ($1):($2) with lines lt 2 title 'Analytical'
pause -1
#set output 'mandeltime_1x10E5.eps'
set title "Pressure Solutions for Mandel Problem \n Time 1.0E5"
plot "dxelast10.ht01/prsr000.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Init Pressr' with lines lt -1 , "dxelast10.ht01/prsr100.stoc" every ::500::599 using ($0/A):(xnorm*$1) title 'Numerical' with linespoints lt 1, "analytical/mandelteor10.5dat" using ($1):($2) with lines lt 2 title 'Analytical'


