set terminal postscript eps "Helvetica" 24
set output 'kk8add-asp-lys.eps'
set grid
set key left top
set xlabel "Debye Length [A]"
set ylabel "Sensitivity"
set xrange[0:9]
set yrange[-0.6:0.6]
plot 'kk8add-asp.dat' title "Asp close to wire" with lp lt 1 lc 1 lw 4 pt 1 ps 2,\
     'kk8add-lys.dat' title "Lys close to wire" with lp lt 1 lc 3 lw 4 pt 2 ps 2
set output
