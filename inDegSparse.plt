set key top right
set logscale xy 10
#set format x "10^{%L}"
set mxtics 10
#set format y "10^{%L}"
set mytics 10
set grid
set xlabel "In-degree"
set ylabel "Count"
set tics scale 2
set xrange [0.75:]
plot 	"inDeg.ModelSparse.tab" using 1:2 title "Model" with linespoints pt 6, "inDeg.KronSmallSparse.tab" using 1:2 title "Kron" with linespoints pt 6
