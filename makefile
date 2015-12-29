#Thesis project makefile
CC=g++
CFLAGS=-c -g -O2 -Wall -Wextra -Isrc -rdynamic -DNDEBUG $(OPTFLAGS)
LIBS = -lgsl -lgslcblas -lglpk -lm $(OPTLIBS)

SRC = $(wildcard src/*.cpp)
OBJ = $(patsubst %.cpp,%.o,$(SRC))

TEST_SRC=$(wildcard tests/*_tests.cpp)
TEST_OBJ=$(patsubst %.cpp,%.o,$(TEST_SRC))
TESTS=$(patsubst %.cpp,%,$(TEST_SRC))

TARGET=build/libparamid.a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

all: main.o $(TARGET) $(SO_TARGET) tests
	$(CC) -g $(OBJ) $< $(LIBS) -o bin/prog

%.o: %.cpp
	$(CC) $(CFLAGS) $< $(LIBS) -o $@

build:
	@mkdir -p build
	@mkdir -p bin
	
dev: CFLAGS=-c -g -Wall -Wextra 
dev: all

$(TARGET): CFLAGS += -fPIC
$(TARGET): build $(OBJ)
	ar rcs $@ $(OBJ)
	
$(SO_TARGET): $(TARGET)
	$(CC) -shared -o $@ $(OBJ)

.PHONY: tests
tests: CFLAGS=-Isrc
tests: LIBS=-Lbuild -lparamid -lgsl -lgslcblas -lboost_system -lboost_unit_test_framework
tests: $(TESTS)
	export LD_LIBRARY_PATH=$(PWD)/build 
#	sh ./tests/runtests.sh

$(TESTS): $(TEST_SRC)
	$(CC) $(CFLAGS) $< $(LIBS) -o $@

clean:
	rm -rf main.o $(OBJ)

cleantests:
	rm -rf $(TESTS)
	
cleanall: clean cleantests
	rm -rf bin build