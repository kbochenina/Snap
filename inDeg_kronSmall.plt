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
set ylabel "Count (non-cum)"
set tics scale 2
plot 	"inDeg.model.tab" using 1:2 title "Model" with linespoints pt 6, "inDeg.kronSmall.tab" using 1:2 title "KronSmall" with linespoints pt 6
