all:
	wave

build: wave.c
	gcc wave.c -o wave -lm

clean:
	rm wave wave_output.dat

run:
	./wave

gnu:
	gnuplot plot_wave.gp