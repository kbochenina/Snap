#set terminal postscript eps enhanced
#set output "test.ps"
set key top right
#set logscale xy 10
#set format x "10^{%L}"
#set mxtics 10
#set format y "10^{%L}"
#set mytics 10
set grid
set xlabel "Clustering coefficient"
set ylabel "Degree"
set tics scale 2
plot 	"ccf.model.tab" using 1:2 title "Model" with linespoints pt 6, "ccf.Small.tab" using 1:2 title "ModelSmall" with linespoints pt 6
