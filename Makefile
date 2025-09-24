all:
	wave

build: wave.c
	gcc wave.c -o wave -lm

clean:
	rm wave timestep/*.txt *.txt

run:
	./wave

gif:
	gnuplot plot_animate.gp