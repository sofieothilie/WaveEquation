set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output 'wave_output.png'

set title "Wave Equation Output"
set xlabel "x"
set ylabel "u(x, t)"
set grid
set key outside

# Set the range for x and y (adjust if needed)
set xrange [0:1]
set yrange [-1:1]

# Animate over time steps
plot 'wave_output.dat' index 0 using 1:2 with lines title "t = 0.004500", \
     '' index 1 using 1:2 with lines title "t = 0.229500", \
     '' index 2 using 1:2 with lines title "t = 0.454500", \
     '' index 3 using 1:2 with lines title "t = 0.679500", \
     '' index 4 using 1:2 with lines title "t = 0.904500", \
     '' index 5 using 1:2 with lines title "t = 1.129500", \
     '' index 6 using 1:2 with lines title "t = 1.354500", \
     '' index 7 using 1:2 with lines title "t = 1.579500", \
     '' index 8 using 1:2 with lines title "t = 1.804500"