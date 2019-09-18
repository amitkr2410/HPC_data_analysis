# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2016 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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

# Include the configuration.
# -include Makefile.inc
# module load root
# source  ${ROOTSYS}/bin/thisroot.sh
PREFIX_BIN=/wsu/home/fy/fy41/fy4125/Software/pythia8230/bin
PREFIX_INCLUDE=/wsu/home/fy/fy41/fy4125/Software/pythia8230/include
PREFIX_LIB=/wsu/home/fy/fy41/fy4125/Software/pythia8230/lib
PREFIX_SHARE=/wsu/home/fy/fy41/fy4125/Software/pythia8230/share/Pythia8
ENABLE_SHARED=true
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so


FASTJET3_USE=true
FASTJET3_BIN=/wsu/home/fy/fy41/fy4125/Software/fastjet/bin/
FASTJET3_INCLUDE=/wsu/home/fy/fy41/fy4125/Software/fastjet/include
FASTJET3_LIB=/wsu/home/fy/fy41/fy4125/Software/fastjet/lib
# Handle GZIP support.
ifeq ($(GZIP_USE),true)
  CXX_COMMON+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
  CXX_COMMON+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
endif

FASTJET3_USE=true
# Check distribution (use local version first, then installed version).
#ifneq ("$(wildcard /wsu/home/fy/fy41/fy4125/Software/pythia8230/lib/libpythia8.a)","")
  PREFIX_LIB=/wsu/home/fy/fy41/fy4125/Software/pythia8230/lib
  PREFIX_INCLUDE=/wsu/home/fy/fy41/fy4125/Software/pythia8230/include
#endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) -I$(PREFIX_INCLUDE)/Pythia8 $(CXX_COMMON) -Wl,-rpath,$(PREFIX_LIB) -ldl

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "Usage: make mainXX"

# The Makefile configuration.
Makefile.inc:
	$(error Error: PYTHIA must be configured, please run "./configure"\
                in the top PYTHIA directory)

# PYTHIA libraries.
$(PREFIX_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

# Examples without external dependencies.
main% : main%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

# GZIP (required).
main34: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(GZIP_USE),true)
	$(CXX) $^ -o $@ $(CXX_COMMON)
else
	@echo "Error: $@ requires GZIP"
endif

# HEPMC2.
main41 main42 main43 main85 main86 main87 main88 main89: $$@.cc\
	$(PREFIX_LIB)/libpythia8.a
ifeq ($(HEPMC2_USE),true)
	$(CXX) $^ -o $@ -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC
else
	@echo "Error: $@ requires HEPMC2"
endif

# PROMC.
main46: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(PROMC_USE),true)
	$(CXX) $^ -o $@ -I$(PROMC_INCLUDE)/src -I$(PROMC_INCLUDE)/include\
	 $(CXX_COMMON) -DPROMC=\"$(PROMC_INCLUDE)\" -Wno-long-long\
	 -L$(PROMC_LIB) -Wl,-rpath,$(PROMC_LIB) -lpromc -lprotoc -lprotobuf\
	 -lprotobuf-lite -lcbook
else
	@echo "Error: $@ requires PROMC"
endif

# EVTGEN (and HEPMC2).
main48: $$@.cc $(PREFIX_LIB)/libpythia8.so
ifeq ($(EVTGEN_USE)$(HEPMC2_USE),truetrue)
	$(CXX) $< -o $@ -I$(EVTGEN_INCLUDE) $(CXX_COMMON)\
	 -DEVTGEN_PYTHIA -DEVTGEN_EXTERNAL -Wl,-rpath,$(HEPMC2_LIB)\
	 -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8\
	 -L$(EVTGEN_LIB) -Wl,-rpath,$(EVTGEN_LIB) -lEvtGenExternal -lEvtGen
else
	@echo "Error: $@ requires EVTGEN and HEPMC2"
endif

# FASTJET3.
main71 main72: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE),true)
	$(CXX) $^ -o $@ -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet
else
	@echo "Error: $@ requires FASTJET3"
endif

# FASTJET3 and HEPMC2.
main81 main82 main83 main84: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE)$(HEPMC2_USE),truetrue)
	$(CXX) $^ -o $@ -I$(FASTJET3_INCLUDE) -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet
else
	@echo "Error: $@ requires FASTJET3 and HEPMC2"
endif

#FastJet with ROOT by Amit
JetpTCrossSectionHadron JetpTCrossSectionHadronRecoil PartonJetpTCrossSectionISR PartonJetpTCrossSection PartonJetpTCrossSectionISRRecoil JetpTCrossSectionPartonHadronRecoil AnalysisSpectra2 AnalysisSpectra3 : $$@.C $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< $(PREFIX_LIB)/libpythia8.a -o $@ -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
         -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet -w -I$(ROOT_INCLUDE) $(CXX_COMMON)\
         `$(ROOTBIN)root-config --cflags` -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`

#Pythia with ROOT
AnalysisSpectraSingleHadron  JetpTCrossSectionPartonHadronPbPb : $$@.C $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< $(PREFIX_LIB)/libpythia8.a -o $@  $(CXX_COMMON)\
          -Wl,-rpath, -w -I$(ROOT_INCLUDE)   `$(ROOTBIN)root-config --cflags` -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`

# ROOT (turn off all warnings for readability).
main91: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(ROOT_USE),true)
	$(CXX) $^ -o $@ -w -I$(ROOT_INCLUDE) $(CXX_COMMON)\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif
main92: $$@.cc $$@.h $$@LinkDef.h $(PREFIX_LIB)/libpythia8.a
ifeq ($(ROOT_USE),true)
	export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT_LIB);\
	 $(ROOT_BIN)rootcint -f $@Dct.cc -c -I$(PREFIX_INCLUDE) $@.h $@LinkDef.h
	$(CXX) $@Dct.cc $^ -o $@ -w -I$(ROOT_INCLUDE) $(CXX_COMMON)\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif

# User-written examples for tutorials, without external dependencies.
mymain% : mymain%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

# Internally used tests, without external dependencies.
test% : test%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

##Adding (Amit)
#HMpair: $$@.cc $(PREFIX_LIB)/libpythia8.a
#       $(CXX) $^ -o $@ $(CXX_COMMON)

#Adding (Amit)
HQpair : HQpair.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)
fgfun	: fgfun.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON) 

fgfun_gluon: fgfun_gluon.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON)

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f weakbosons.lhe; rm -f Pythia8.promc; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*
