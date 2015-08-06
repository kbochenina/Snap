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
#set terminal png size 1000,800
#set output 'inDeg.Model.png'
plot 	"outDegC.modelSmall.tab" using 1:2 title "Model" with linespoints pt 6, "outDegC.scale0.tab" using 1:2 title "scale0" with linespoints pt 6, "outDegC.scale1.tab" using 1:2 title "scale1" with linespoints pt 6, "outDegC.scale2.tab" using 1:2 title "scale2" with linespoints pt 6, "outDegC.scale3.tab" using 1:2 title "scale3" with linespoints pt 6, "outDegC.scale4.tab" using 1:2 title "scale4" with linespoints pt 6, "outDegC.scale5.tab" using 1:2 title "scale5" with linespoints pt 6, "outDegC.scale6.tab" using 1:2 title "scale6" with linespoints pt 6, "outDegC.scale7.tab" using 1:2 title "scale7" with linespoints pt 6, "outDegC.scale8.tab" using 1:2 title "scale8" with linespoints pt 6, "outDegC.scale9.tab" using 1:2 title "scale9" with linespoints pt 6
