set key bottom right
set grid
set logscale xy
set xlabel "Number of hops"
set ylabel "Reachable pairs of nodes"
plot 	"hop.model.tab" using 1:2 title "Model" with linespoints pt 6, "hop.kronSingle.tab" using 1:2 title "KronSingle" with linespoints pt 6
