set key bottom right
set grid
set logscale xy
set xlabel "Number of hops"
set ylabel "Reachable pairs of nodes"
plot 	"hop.modelsmall.tab" using 1:2 title "Model" with linespoints pt 6, "hop.kronsmall.tab" using 1:2 title "Kron" with linespoints pt 6
