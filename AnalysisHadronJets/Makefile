# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2023 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash
#source /snap/root-framework/954/usr/local/bin/thisroot.sh
# Include the configuration.
#-include Makefile.inc
PYTHIA=/nfs/zfs2/Software/pythia8309
PREFIX_BIN=${PYTHIA}/bin
PREFIX_INCLUDE=${PYTHIA}/include
PREFIX_LIB=${PYTHIA}/lib
PREFIX_SHARE=${PYTHIA}/share/Pythia8
ENABLE_SHARED=true
CXX=g++  -std=c++17
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so


CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl


#ROOT_LIB=/snap/root-framework/954/usr/local/lib
#ROOT_INCLUDE=/snap/root-framework/954/usr/local/include
#ROOT_BIN=/snap/root-framework/954/usr/local/include/bin




################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	$(info Usage: make mainXX)

# PYTHIA library.
ROOT_INCLUDE=/snap/root-framework/950/usr/local/include
ROOT_LIB=/snap/root-framework/950/usr/local/lib
ROOT_CFLAGS:=$(shell root-config --cflags)
ROOT_LIBS:=$(shell root-config --libs)

#ROOT
% : %.C
	$(CXX) $< $(PREFIX_LIB)/libpythia8.a -o $@  $(CXX_COMMON) \
	-Wl,-rpath, -w  $(ROOT_CFLAGS) -I$(PREFIX_INCLUDE) \
	-L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl $(ROOT_LIBS)

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat;\
	rm -f weakbosons.lhe; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;\
	rm -f *.log; rm -f *plot.py; rm -f *.pcm; rm -f *.mpi;\
	rm -f *.hepmc; rm -f *.yoda; rm -f *.root;
