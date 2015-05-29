set key top right
set logscale xy 10
#set format x "10^{%L}"
set mxtics 10
#set format y "10^{%L}"
set mytics 10
set grid
set xlabel "Rank"
set ylabel "Eigenvalue"
set tics scale 2
set style line 1 lt 1 lc 1 lw 1 pt 1
set style line 2 lt 1 lc 7 lw 1 pt 1
set style line 3 lt 1 lc 3 lw 1 pt 1
set style line 4 lt 2 lc 1 lw 1 pt 4
set style line 5 lt 2 lc 7 lw 1 pt 4
set style line 6 lt 2 lc 3 lw 1 pt 4
#set terminal png size 1000,800
#set output 'inDeg.Model.png'
plot 	"eigval.ModelEigen0.tab" using 1:2 title "ModelEigen0" with linespoints ls 1, "eigval.ModelEigen1.tab" using 1:2 title "ModelEigen1" with linespoints ls 2, "eigval.ModelEigen2.tab" using 1:2 title "ModelEigen2" with linespoints ls 3, "ApproxEigen0.tab" using 1:2 title "ApproxEigen0" with linespoints ls 4, "ApproxEigen1.tab" using 1:2 title "ApproxEigen1" with linespoints ls 5, "ApproxEigen2.tab" using 1:2 title "ApproxEigen2" with linespoints ls 6
