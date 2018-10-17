set terminal png
set output 'pic.png'
plot 'fort.77' u 1:2 w l
set xlabel [0,4]
