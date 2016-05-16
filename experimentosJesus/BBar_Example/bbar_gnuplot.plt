##
##  FOR GNUPLOT VERSION 5.XX uncomment next line
#set colorsequence classic
##
reset

#
## Define Geometric and Inertia Values
#
P=1.0
E=1.0
c=2.0
L=16.0
I=2.0*(c**3)/3.0
#
## Define configurations for graphics display
#
set grid
#set zeroaxis lt -1 lw 3
set pointsize 2
set xlabel "X Axis"
set ylabel "u_y(x,0)"
set key bottom left
set size square
set xrange [-0:16]
#set yrange [-400:1.0]
#
## Define output parameters
#set terminal postscript eps enhanced color solid
#monochrome dashed
#
## Setup Poissson Value
#
poisson=0.300
#
KK=-P*(1.0-poisson**2)/(6.0*E*I)
KKlin=(4.0+5.0*(poisson/(1.0-poisson)))
#
#set output 'nu300c.eps'
#
## Plot de functions
#
## Set Title
#
set title 'Solutions for Vertical Displacements. Poisson=0.3'
plot "bbar_nu.300/displac.dat" every 5 using ($0):($2) title 'B-Bar Formulation' with linespoints lt 1 lw 1, KK*(KKlin*(c**2)*x+(3.0*L-x)*x**2) title 'Analytic Solution' with lines lt -1 lw 1, "classic_nu.300/displac.dat" every 5 using ($0):($2) title 'Classic Formulation' with linespoints lt 2 lw 1
pause -1
#
## Setup Poissson Value
#
poisson=0.450
#
KK=-P*(1.0-poisson**2)/(6.0*E*I)
KKlin=(4.0+5.0*(poisson/(1.0-poisson)))
#
#set output 'nu450c.eps'
#
## Plot de functions
#
## Set Title
#
set title 'Solutions for Vertical Displacements. Poisson =0.45'
plot "bbar_nu.450/displac.dat" every 5 using ($0):($2) title 'B-Bar Formulation' with linespoints lt 1 lw 1, KK*(KKlin*(c**2)*x+(3.0*L-x)*x**2) title 'Analytic Solution' with lines lt -1 lw 1,  "classic_nu.450/displac.dat" every 5 using ($0):($2) title 'Classic Formulation' with linespoints lt 2 lw 1
pause -1
## Setup Poissson Value
#
poisson=0.475
#
KK=-P*(1.0-poisson**2)/(6.0*E*I)
KKlin=(4.0+5.0*(poisson/(1.0-poisson)))
#
#
#set output 'nu475c.eps'
#
## Plot de functions
#
## Set Title
#
set title 'Solutions for Vertical Displacements. Poisson =0.475'
plot "bbar_nu.475/displac.dat" every 5 using ($0):($2) title 'B-Bar Formulation' with linespoints lt 1 lw 1, KK*(KKlin*(c**2)*x+(3.0*L-x)*x**2) title 'Analytic Solution' with lines lt -1 lw 1,  "classic_nu.475/displac.dat" every 5 using ($0):($2) title 'Classic Formulation' with linespoints lt 2 lw 1
#pause -1
#
## Setup Poissson Value
#
poisson=0.485
#
KK=-P*(1.0-poisson**2)/(6.0*E*I)
KKlin=(4.0+5.0*(poisson/(1.0-poisson)))
#
#
#set output 'nu485c.eps'
#
## Plot de functions
#
## Set Title
#
set title 'Solutions for Vertical Displacements. Poisson =0.485'
plot "bbar_nu.485/displac.dat" every 5 using ($0):($2) title 'B-Bar Formulation' with linespoints lt 1 lw 1, KK*(KKlin*(c**2)*x+(3.0*L-x)*x**2) title 'Analytic Solution' with lines lt -1 lw 1,  "classic_nu.485/displac.dat" every 5 using ($0):($2) title 'Classic Formulation' with linespoints lt 2 lw 1
#pause -1
#
#set output 'nu490c.eps'
#
## Setup Poissson Value
#
poisson=0.490
#
KK=-P*(1.0-poisson**2)/(6.0*E*I)
KKlin=(4.0+5.0*(poisson/(1.0-poisson)))
#
## Plot de functions
#
## Set Title
#
set title 'Solutions for Vertical Displacements. Poisson =0.490'
plot "bbar_nu.490/displac.dat" every 5 using ($0):($2) title 'B-Bar Formulation' with linespoints lt 1 lw 1, KK*(KKlin*(c**2)*x+(3.0*L-x)*x**2) title 'Analytic Solution' with lines lt -1 lw 1,  "classic_nu.490/displac.dat" every 5 using ($0):($2) title 'Classic Formulation' with linespoints lt 2 lw 1
pause -1
#
#set output 'nu499c.eps'
#
## Setup Poissson Value
#
#poisson=0.499
#
KK=-P*(1.0-poisson**2)/(6.0*E*I)
KKlin=(4.0+5.0*(poisson/(1.0-poisson)))
#
## Plot de functions
#
## Set Title
#
set title 'Solutions for Vertical Displacements. Poisson =0.499'
plot "bbar_nu.499/displac.dat" every 5 using ($0):($2) title 'B-Bar Formulation' with linespoints lt 1 lw 1, KK*(KKlin*(c**2)*x+(3.0*L-x)*x**2) title 'Analytic Solution' with lines lt -1 lw 1,  "classic_nu.499/displac.dat" every 5 using ($0):($2) title 'Classic Formulation' with linespoints lt 2 lw 1

