all:
	wave

buildseq: wave.c
	gcc wave.c -o wave -lm

buildpar:
	nvcc -O3 -o wave wave_parallelised.cu

clean:
	rm wave timestep/*.txt ; rmdir timestep

run:
	./wave

gif:
	gnuplot plot_animate.gp

seq: clean buildseq run gif

par: clean buildpar run gif