set terminal gif animate delay 10
set output 'wave3d.gif'

set view 60, 30
set xrange [-1:1]
set yrange [-1:1]
set zrange [-1:1]
set cbrange [-1:1]
set palette defined (0 "blue", 1 "white", 2 "red")

do for [i=0:19] {
    fname = sprintf("timestep/wave3d_%03d.txt", i)
    eval sprintf("splot '%s' using 1:2:3:4 with points palette pointsize 0.3 notitle", fname)
}
