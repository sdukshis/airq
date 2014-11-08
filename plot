#!/usr/bin/env gnuplot -persist
#set terminal postscript eps enhanced  dashed font "/usr/share/fonts/truetype/times.ttf" 14
set encoding koi8r
set output "ALL.eps"
set terminal postscript enhanced "Courier" 14
set size ratio 0.5
set xlabel "r_{min}" 
set ylabel "C/C_{AAR}"
set y2label "" 
set yrange [0.8:1]
set xrange [8000:13000]
set xtics 8000,1000,14000
set ytics nomirror
set ytics 0.8, .05, 1
set mytics 2
#set y2tics 5
set my2tics 2
set grid
#set autoscale y
set key bottom right
set key box linestyle 1
#plot "results/SVO.res" using 1:2 smooth csplines title "SVO"  with lines lt -1
#plot "results/DME.res" using 1:2 smooth csplines title "DME"  with lines lt -1
plot "results/SVO1_1.res" title "SVO07" with linespoints, \
     "results/SVO1_2.res" title "SVO25" with linespoints, \
     "results/DME1_1.res" title "DME09" with linespoints, \
     "results/DME1_2.res" title "DME32" with linespoints, \
     "results/VKO1_1.res" title "VKO01" with linespoints, \
     "results/VKO1_2.res" title "VKO19" with linespoints
