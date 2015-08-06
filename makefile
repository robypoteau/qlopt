# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall -g -Wextra -DNDEBUG
OBJ = main.o thesis_functions.o numerical_integration.o nonlinear_odes.o spline.o bspline.o
LIBS = -lmpfr -lgmp -lgsl -lgslcblas

all: prog

prog: $(OBJ)
	$(CC) -g $(OBJ) $(LIBS) -o prog
	
main.o: main.cpp
	$(CC) $(CFLAGS) $(LIBS) main.cpp

thesis_functions.o: thesis_functions.cpp
	$(CC) $(CFLAGS) $(LIBS) thesis_functions.cpp

numerical_integration.o: numerical_integration.cpp
	$(CC) $(CFLAGS) $(LIBS) numerical_integration.cpp

nonlinear_odes.o: nonlinear_odes.cpp
	$(CC) $(CFLAGS) $(LIBS) nonlinear_odes.cpp

spline.o: spline.cpp
	$(CC) $(CFLAGS) $(LIBS) spline.cpp
	
bspline.o: bspline.cpp
	$(CC) $(CFLAGS) $(LIBS) bspline.cpp	

dev: CFLAGS=-c -g -Wall -Wextra
dev: all

clean:
	rm -rf *o prog