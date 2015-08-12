#Thesis project makefile
CC=g++
CFLAGS=-g -O2 -Wall -Wextra -Isrc -rdynamic -DNDEBUG $(OPTFLAGS)
LIBS = -lmpfr -lgmp -lgsl -lgslcblas $(OPTLIBS)

SRC = $(wildcard *.c)
OBJ = main.o thesis_functions.o numerical_integration.o nonlinear_odes.o spline.o bspline.o latex_output.o

all: prog

test:
	@echo $(SRC)

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

latex_output.o: latex_output.c
	$(CC) $(CFLAGS) $(LIBS) latex_output.c	
	
dev: CFLAGS=-c -g -Wall -Wextra
dev: all

clean:
	rm -rf *o prog

