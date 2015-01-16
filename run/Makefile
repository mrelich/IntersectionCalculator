SHELL := /bin/bash

RCXX=$(shell root-config --cflags)
RLIB=$(shell root-config --libs)

CXX		= g++
CFLAGS		= -pipe -g \
-D __SKIP_CLASS_DEF__ \
-Wall -Wextra -Winit-self -Wcast-align -Wno-comment -Wno-ignored-qualifiers \
-O2 -march=native -mfpmath=sse -fPIC \
-Wl,--no-as-needed

#---- Paths ----#
inc = ../include
src = ../src

#---- Links ----#
links = -I$(inc) -I$(ROOTSYS)/include

CXXFLAGS	= $(CFLAGS) $(RCXX)
LDFLAGS		= -lm -L. $(RLIB)

#---- Targets ----#
Targets = test

#---- Generic Commands ----#
all: Makefile $(Targets)

obj/%.o:
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(links)

#---- Additional Objects ----#
obj/Calculator.o:	$(src)/Calculator.cxx
obj/main.o:		$(src)/main.cxx

test:				\
	obj/Calculator.o	\
	obj/main.o
	$(CXX) $+ -o $@ $(CXXFLAGS) $(LDFLAGS) $(links)

clean:
	rm -f $(Targets)
	rm -f obj/*.o
	rm -f obj/*.so
	rm -f obj/*.cxx
	rm -f obj/*.h
