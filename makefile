# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall -g -lmpfr -lgmp

all: prog

prog: main.o thesis_functions.o numerical_integration.o nonlinear_odes.o regs.o
	$(CC) -g main.o thesis_functions.o numerical_integration.o nonlinear_odes.o regs.o -lmpfr -lgmp -o prog
	
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

thesis_functions.o: thesis_functions.cpp
	$(CC) $(CFLAGS) thesis_functions.cpp

numerical_integration.o: numerical_integration.cpp
	$(CC) $(CFLAGS) numerical_integration.cpp

nonlinear_odes.o: nonlinear_odes.cpp
	$(CC) $(CFLAGS) nonlinear_odes.cpp

regs.o: regs.cpp
	$(CC) $(CFLAGS) regs.cpp
	
clean:
	rm -rf *o prog