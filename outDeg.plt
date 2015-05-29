set key bottom right
set logscale xy 10
#set format x "10^{%L}"
set mxtics 10
#set format y "10^{%L}"
set mytics 10
set grid
set xlabel "Out-degree"
set ylabel "Count"
set tics scale 2
plot 	"outDeg.Model.tab" using 1:2 title "" with linespoints pt 6, "outDeg.KronSmall.tab" using 1:2 title "" with linespoints pt 6
