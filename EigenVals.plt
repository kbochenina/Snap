set key bottom right
#set logscale xy 10
#set format x "10^{%L}"
set mxtics 10
#set format y "10^{%L}"
set mytics 10
set grid
set xlabel "Rank"
set ylabel "Eigenvalue"
set tics scale 2
#set terminal png size 1000,800
#set output 'inDeg.Model.png'
plot 	"eigval.ModelEigen.tab" using 1:2 title "Model" with linespoints pt 6, "eigval.KronEigen.tab" using 1:2 title "Kron" with linespoints pt 6, "eigval.ProbMtx.tab" using 1:2 title "Kron" with linespoints pt 6
