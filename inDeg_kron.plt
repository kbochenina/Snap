#set terminal postscript eps enhanced
#set output "test.ps"
set key top right
set logscale xy 10
set mxtics 10
set mytics 10
set grid
set xlabel "In-degree"
set ylabel "Count (non-cum)"
set tics scale 2
plot 	"inDeg.modelmodel.tab" using 1:2 title "Model" with linespoints pt 6, "inDeg.kronmodel.tab" using 1:2 title "Kron" with linespoints pt 6, "inDeg.kron.tab" using 1:2 title "Kron" with linespoints pt 6, "inDeg.kronSmall2.tab" using 1:2 title "Kron" with linespoints pt 6