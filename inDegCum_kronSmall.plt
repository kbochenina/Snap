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
plot 	"inDegC.modelSmall.tab" using 1:2 title "Model" with linespoints pt 6, "inDegC.KronSmall.tab" using 1:2 title "Kron" with linespoints pt 6, "inDegC.KronSmall1.tab" using 1:2 title "Kron" with linespoints pt 6, "inDegC.KronSmall2.tab" using 1:2 title "Kron" with linespoints pt 6
