PWD=$(shell pwd)
UNAME=$(shell uname -s)

# Project name
PRJNAME=paramid

# Variables related to directory structure
BUILDDIR?=$(PWD)/build
BINDIR=$(PWD)/bin
SRCDIR=$(PWD)/src
APPSRCDIR=$(PWD)/examples
# Extension for source files
SRCEXT=cpp

# Extension for dynamic libraries depends on the platform
ifeq ($(UNAME), Linux)
DYLIBEXT:=so
endif
ifeq ($(UNAME), Darwin)
DYLIBEXT:=dylib
endif

# Build configuration
CC=g++
CXXFLAGS_EXTRA?=
CXXFLAGS=-g -O2 -Wall -Wextra -I$(PWD)/src -rdynamic -fPIC -DNDEBUG -std=c++11 $(CXXFLAGS_EXTRA)
LIBS=-L$(BUILDDIR) -l$(PRJNAME)
#LIBS=-lgsl -lgslcblas -lmpfr -lgmp -lm

#TODO: Document variables

APPSRC=$(wildcard $(APPSRCDIR)/*.$(SRCEXT)) # Find all source files under APPSRCDIR
EXE=$(patsubst $(APPSRCDIR)/%.$(SRCEXT),$(BINDIR)/%,$(APPSRC))

SRC=$(wildcard $(SRCDIR)/*.$(SRCEXT)) # Find all source files under SRCDIR
OBJ=$(patsubst $(SRCDIR)/%.$(SRCEXT),$(BUILDDIR)/%.o,$(SRC))

TEST_SRC=$(wildcard tests/*_tests.cpp) # Find all source files under tests/
TEST_OBJ=$(patsubst %.$(SRCEXT),%.o,$(TEST_SRC))
TESTS=$(patsubst %.$(SRCEXT),%,$(TEST_SRC))

TARGET=$(BUILDDIR)/lib$(PRJNAME).a
SO_TARGET=$(patsubst %.a,%.$(DYLIBEXT),$(TARGET))

usage:
	@echo "usage: make [TARGET]"
	@echo "    where TARGET is one of"
	@echo "    - usage (show this message)"
	@echo "    - static: Generate the QLopt static library"
	@echo "    - shared: Generate the QLopt shared library"
	@echo "    - apps: Generate the QLopt executables (under $(APPSRCDIR))"
	@echo "    - tests: Run the included test script (disabled)"

all: static shared apps

#TODO: use common flags for different builds/targets, or document why they need
# to vary.
dev: CXXFLAGS+=-c
dev: all

apps: $(EXE)

# The rule below uses a _pattern_: objects matching the string on the left-hand
# side of the colon (after expanding any variables) using "%" as a wildcard will
# be built from the corresponding object on the right-hand side using the shell
# commands below the rule.
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) -c $(CXXFLAGS) $<  -o $@

$(BINDIR)/%: $(APPSRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CXXFLAGS) $< -Wl,-rpath,$(BUILDDIR) $(LIBS) -o $@

#TODO convert this into a pattern (?) rule also
$(TESTS): $(TEST_SRC)
	$(CC) $(CXXFLAGS) $< $(LIBS) -o $@

.PHONY: static
static: $(TARGET)
$(TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	ar rcs $@ $^

.PHONY: shared
shared: $(SO_TARGET)
$(SO_TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	$(CC) -shared -o $@ $^

.PHONY: tests
tests: CXXFLAGS=-I$(BUILDDIR) -I$(SRCDIR)
tests: $(TESTS)
	export LD_LIBRARY_PATH=$(PWD)/$(BUILDDIR)
#	sh ./tests/runtests.sh

.PHONY: clean
clean:
	rm -rf $(BINDIR) $(BUILDDIR) $(TESTS)
