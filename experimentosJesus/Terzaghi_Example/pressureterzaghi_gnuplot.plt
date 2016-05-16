reset
set grid
set yrange [0:12]
set xrange [0:10.5]
set key top horizontal
set xlabel "Vertical Axis"
set ylabel "Pressure"
set title "Numerical and Analytical Solutions for Time =0.001" 
plot "dxelast10.ht01/prsr001.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step001.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.002" 
plot "dxelast10.ht01/prsr002.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step002.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.003" 
plot "dxelast10.ht01/prsr003.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step003.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.004" 
plot "dxelast10.ht01/prsr004.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step004.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.005" 
plot "dxelast10.ht01/prsr005.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step005.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.006" 
plot "dxelast10.ht01/prsr006.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step006.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.007" 
plot "dxelast10.ht01/prsr007.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step007.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.008" 
plot "dxelast10.ht01/prsr008.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step008.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.009" 
plot "dxelast10.ht01/prsr009.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step009.dat" using ($2):($3) with lines title 'theoric'
pause -1
set title "Numerical and Analytical Solutions for Time =0.010" 
plot "dxelast10.ht01/prsr010.stoc" every 100::50::5000 using (($0+1)/5):($1) with lp title 'numeric', "analytical/prsr_step010.dat" using ($2):($3) with lines title 'theoric'

