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

# Compilation configuration
CXXFLAGS_EXTRA?=
CXXFLAGS=-g -O2 -Wall -Wextra -I$(PWD)/src -I/usr/local/include -rdynamic -fPIC -DNDEBUG -std=c++11 $(CXXFLAGS_EXTRA)

# Linking configuration
LIBS=-L$(BUILDDIR) -l$(PRJNAME) -lsundials_sunnonlinsolnewton -lsundials_sunnonlinsolfixedpoint -lsundials_cvodes -lsundials_nvecserial -L/usr/local/lib

#TODO: Document variables

# Build with `make apps`
APPSRC=$(wildcard $(APPSRCDIR)/*.$(SRCEXT)) # Find all source files under APPSRCDIR
EXE=$(patsubst $(APPSRCDIR)/%.$(SRCEXT),$(BINDIR)/%,$(APPSRC))

# SRC files are compiled into corresponding OBJ objects
SRC=$(wildcard $(SRCDIR)/*.$(SRCEXT)) # Find all source files under SRCDIR
OBJ=$(patsubst $(SRCDIR)/%.$(SRCEXT),$(BUILDDIR)/%.o,$(SRC))

# The variables below are not currently used.
TEST_SRC=$(wildcard tests/*_tests.cpp) # Find all source files under tests/
TEST_OBJ=$(patsubst %.$(SRCEXT),%.o,$(TEST_SRC))
TESTS=$(patsubst %.$(SRCEXT),%,$(TEST_SRC))

# Static library and corresponding shared library
STATIC_TARGET=$(BUILDDIR)/lib$(PRJNAME).a
SHARED_TARGET=$(patsubst %.a,%.$(DYLIBEXT),$(TARGET))

# Phony targets don't result in a file appearing in the filesystem
# This one prints a usage message if you type `make` without any
# argument, or if you type `make usage`.

.PHONY: usage
usage:
	@echo "usage: make [TARGET]"
	@echo "    where TARGET is one of"
	@echo "    - usage (show this message)"
	@echo "    - static: Generate the QLopt static library"
	@echo "    - shared: Generate the QLopt shared library"
	@echo "    - apps: Generate the QLopt executables (from source files in $(APPSRCDIR))"
	@echo "    - tests: Run the included test script (disabled)"

.PHONY: all
all: static shared apps

.PHONY: apps
apps: $(EXE)

# The rule below uses a _pattern_: objects matching the string on the left-hand
# side of the colon (after expanding any variables) using "%" as a wildcard will
# be built from the corresponding object on the right-hand side using the shell
# commands below the rule.
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)

	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BINDIR)/%: $(APPSRCDIR)/%.$(SRCEXT) $(SHARED_TARGET)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $< -Wl,-rpath,$(BUILDDIR) $(LIBS) -o $@

.PHONY: static
static: $(STATIC_TARGET)
$(STATIC_TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	ar rcs $@ $^


.PHONY: shared
shared: $(SHARED_TARGET)
$(SHARED_TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	$(CXX) -shared -o $@ $^

.PHONY: tests
tests: $(EXE)
	$(foreach single_exe,$(EXE),$(single_exe);)

.PHONY: clean
clean:
	rm -rf $(BINDIR) $(BUILDDIR) $(TESTS)
