set key bottom right
set logscale xy 10
#set format x "10^{%L}"
set mxtics 10
#set format y "10^{%L}"
set mytics 10
set grid
set xlabel "Clustering coefficient"
set ylabel "Node degree"
set tics scale 2
#set terminal png size 1000,800
#set output 'inDeg.Model.png'
plot 	"ccf.Model.tab" using 1:2 title "Model" with linespoints pt 6, "ccf.KronSingle.tab" using 1:2 title "Kron" with linespoints pt 6
