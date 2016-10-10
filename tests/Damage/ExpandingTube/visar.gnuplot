plot "VISAR-a" using   (1e6*$2 - 14.5):(0.01*$3) title "Probe A" with linespoints lt 0 ps 2
replot "VISAR-b" using (1e6*$2 - 14.5):(0.01*$3) title "Probe B" with linespoints lt 1 ps 2
replot "VISAR-c" using (1e6*$2 - 14.5):(0.01*$3) title "Probe C" with linespoints lt 2 ps 2

replot "VISAR-experimental/25_mm.txt" using (1e6*$1):(1e3*$2) title "Probe A -- Experiment" with lines lt 0
replot "VISAR-experimental/20_mm.txt" using (1e6*$1):(1e3*$2) title "Probe B -- Experiment" with lines lt 1
replot "VISAR-experimental/15_mm.txt" using (1e6*$1):(1e3*$2) title "Probe C -- Experiment" with lines lt 2

set xrange [0:15]
set yrange [-120:]

set key bottom right
#set title "VISAR velocity probes"
set xlabel "time (microsec)"
set ylabel "velocity (meters/sec)"
