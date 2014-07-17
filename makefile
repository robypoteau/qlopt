# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall -g -lmpfr -lgmp

all: prog

prog: main.o nonlinear_odes.o misc.o numerical_integration.o thesis_functions.o
	$(CC) -g main.o nonlinear_odes.o misc.o numerical_integration.o thesis_functions.o -lmpfr -lgmp -o prog
	
main.o: main.cpp
	$(CC) $(CFLAGS)  main.cpp

nonlinear_odes.o: nonlinear_odes.cpp
	$(CC) $(CFLAGS) -I /mnt/hgfs/thesis/src/ nonlinear_odes.cpp

misc.o: misc.cpp
	$(CC) $(CFLAGS) misc.cpp
	
numerical_integration.o: numerical_integration.cpp
	$(CC) $(CFLAGS) numerical_integration.cpp

thesis_functions.o: thesis_functions.cpp
	$(CC) $(CFLAGS) thesis_functions.cpp
	
clean:
	rm -rf *o prog