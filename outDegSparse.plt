set key top right
set logscale xy 10
#set format x "10^{%L}"
set mxtics 10
#set format y "10^{%L}"
set mytics 10
set grid
set xlabel "Out-degree"
set ylabel "Count"
set tics scale 2
plot 	"outDeg.ModelSparse.tab" using 1:2 title "Model" with linespoints pt 6, "outDeg.KronModelSparse.tab" using 1:2 title "KronModel" with linespoints pt 6
