#!/bin/sh

awk '/Size/ { print $2 " " $4 }' timing-$1.out > timing-$1.dat

cat > timing.gp <<EOF
set term postscript
set output 'timing-$1.eps'
set xlabel 'Dimension'
set ylabel 'Mflop/s'
set nokey
set grid 
plot 'timing-$1.dat' using 1:2 with lines
EOF

gnuplot timing.gp
rm -f timing.gp
epstopdf timing-$1.eps
