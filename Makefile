all:
	wave

buildseq: wave.c
	gcc wave.c -o wave -lm

buildpar:
	nvcc -O3 -o wave wave_parallelised.cu

clean:
	rm wave timestep/*.txt

run:
	./wave

gif:
	gnuplot plot_animate.gp