#
# This is a makefile for lab 1.
#

CC = g++
LINK_OPTION = -lz
COMPILE_OPTION = -g

all: sim

clean:
	rm -f *.o *.out sim

sim: all_knobs.o knob.o
	${CC} -o sim sim.cpp all_knobs.o knob.o ${LINK_OPTION}

all_knobs.o: all_knobs.cpp all_knobs.h
	${CC} -c ${COMPILE_OPTION} all_knobs.cpp

knob.o: knob.cpp knob.h
	${CC} -c ${COMPILE_OPTION} knob.cpp
