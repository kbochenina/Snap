#set terminal postscript eps enhanced
#set output "test.ps"
set key bottom right
#set logscale y 10
#set format x "10^{%L}"
#set mxtics 10
#set format y "10^{%L}"
#set mytics 10
set grid
set xlabel "Number of hops"
set ylabel "Reachable pairs of nodes"
# set tics scale 2
plot 	"hop.model.tab" using 1:2 title "Model" with linespoints pt 6, "hop.kronmodel.tab" using 1:2 title "KronModel" with linespoints pt 6