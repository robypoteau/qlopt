#Dissertation research makefile
CC=g++
CXXFLAGS=-g -O2 -Wall -Wextra -Isrc -Lbuild -rdynamic -DNDEBUG -std=c++11
LIBS=-lgsl -lgslcblas -lmpfr -lgmp -lm
PRJNAME=paramid

BUILDDIR=build
BINDIR=bin
SRCDIR=src
APPSRCDIR=appsrc
SRCEXT=cpp

APPSRC=$(shell find $(APPSRCDIR) -type f -name "*.$(SRCEXT)")
EXE=$(patsubst $(APPSRCDIR)/%,$(BINDIR)/%,$(APPSRC:.$(SRCEXT)=))

SRC=$(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJ=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SRC:.$(SRCEXT)=.o))

TEST_SRC=$(wildcard tests/*_tests.cpp)
TEST_OBJ=$(patsubst %.cpp,%.o,$(TEST_SRC))
TESTS=$(patsubst %.cpp,%,$(TEST_SRC))

TARGET=$(BUILDDIR)/lib$(PRJNAME).a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

all: $(TARGET) $(SO_TARGET) $(EXE) tests

dev: CXXFLAGS=-c -g -Wall -Wextra -std=c++11
dev: all

apps: $(BINDIR) $(EXE)

#Redefining rules, this with the % refine rules
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CC) -c $(CXXFLAGS) $< $(LIBS) -o $@

$(BINDIR)/%:$(APPSRCDIR)/%.$(SRCEXT)
	$(CC) $(CXXFLAGS) $< -l$(PRJNAME) -Wl,-rpath=build $(LIBS) -o $@
#End of redfinitions

#TODO convert this into a refefinition rule also
$(TESTS): $(TEST_SRC)
	$(CC) $(CXXFLAGS) $< $(LIBS) -o $@

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

$(BINDIR):
	@mkdir -p $(BINDIR)

$(TARGET): CXXFLAGS += -fPIC
$(TARGET): $(BUILDDIR) $(BINDIR) $(OBJ)
	ar rcs $@ $(OBJ)

$(SO_TARGET): $(TARGET)
	$(CC) -shared -o $@ $(OBJ)

.PHONY: tests
tests: CXXFLAGS=Ibuild
tests: LIBS=-L$(BUILDDIR) -l$(PRJNAME) -lgsl -lgslcblas -lboost_system -lboost_unit_test_framework
tests: $(TESTS)
	export LD_LIBRARY_PATH=$(PWD)/$(BUILDDIR)
#	sh ./tests/runtests.sh


cleanapps:
	rm -rf $(BINDIR)

cleanbuild:
	rm -rf $(BUILDDIR)

cleantests:
	rm -rf $(TESTS)
