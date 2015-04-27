#set terminal postscript eps enhanced
#set output "test.ps"
set key top right
set logscale xy 10
set format x "10^{%L}"
set mxtics 10
set format y "10^{%L}"
set mytics 10
set grid
set xlabel "In-degree"
set ylabel "Count (cum)"
set tics scale 2
plot 	"inDegC.modelSparse.tab" using 1:2 title "Model" with linespoints pt 6, "inDegC.kronmodelSparse.tab" using 1:2 title "KronModel" with linespoints pt 6
