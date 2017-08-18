############################################################
# Makefile for DDCalc                                      #
#                                                          #
# Simple usage:                                            #
#   make all                                               #
#                                                          #
#   GAMBIT Dark Matter Workgroup                           #
#   ddcalc@projects.hepforge.org                           #
#                                                          #
############################################################

#################### INITIALIZATION ########################

# Generate variables needed for substitution.
empty :=
space := ${empty} ${empty}
lpar  := ${empty}(${empty}
rpar  := ${empty})${empty}

#################### VERSION ###############################

# Software release version.  Extracted from main code file.
DDDCALC_VERSION=$(shell grep "VERSION_STRING =" $(SRC)/DDConstants.f90 | cut -d \' -f2)
#$(info DDCalc version is $(DDDCALC_VERSION))

#################### DIRECTORIES ###########################

DDCALC_DIR := $(shell pwd)
SRC := $(DDCALC_DIR)/src
INCLUDE := $(DDCALC_DIR)/include
EXAMPLES := $(DDCALC_DIR)/examples
DOC := $(DDCALC_DIR)/doc
DATA := $(DDCALC_DIR)/data
BUILD := $(DDCALC_DIR)/build
LIB := $(DDCALC_DIR)/lib

#################### COMPILER ##############################

# The fortran compiler is set here, along with the
# compilation options.  If the compiler is not specified via
# an environmental variable or as an argument to the make
# invocation (FC=<...>), make will search for ifort and
# gfortran, in that order, and use one of those.

# If you want to set your own FFLAGS, set FOPT and this
# makefile will ensure that the appropriate things are appended.

# Alternate compilation flags intended for debugging can be
# used by setting DEBUG to any non-empty value, which may be
# done from the command line with e.g. 'make DEBUG=Y'.

# If FC is not set, search for ifort or gfortran compilers,
# in that order.  If neither is found, make stops with an
# error message.
ifneq (,$(findstring "$(origin FC)","undefined" "default"))
  FC := $(or $(shell which ifort 2> /dev/null),\
             $(shell which gfortran 2> /dev/null),\
             $(error Could not find ifort or gfortran compiler \
                  in the path.  A valid compiler must be given.))
  FC := $(notdir $(FC))
endif

# Check if using ifort or gfortran and give a warning
# otherwise.  These are the only supported compilers.
ifeq (,$(findstring ifort,$(FC))$(findstring gfortran,$(FC)))
  $(info $(empty)------------------------------------------------------------)
  $(info $(empty)  WARNING:)
  $(info $(empty)  Only the ifort $(lpar)Intel$(rpar) and gfortran \
                 $(lpar)GCC$(rpar) compilers have)
  $(info $(empty)  been tested with DDCalc.  Proceeding with unsupported)
  $(info $(empty)  compiler $(FC).)
  $(info $(empty)------------------------------------------------------------)
  #$(info $(empty))
endif

# Fortran compiler flags.  Update only if not already specified.
ifneq (,$(findstring "$(origin FOPT)","undefined" "default"))
  # Default flags.
  ifeq (,$(DEBUG))
    ifneq (,$(findstring ifort,$(FC)))
      FFLAGS := -fast -module $(BUILD) -fpp
    else ifneq (,$(findstring gfortran,$(FC)))
      FFLAGS := -O3 -fno-range-check -J $(BUILD) -cpp
    else
      $(info WARNING: No default compilation flags for compiler $(FC).)
    endif
  # Debugging flags, activated by setting DEBUG non-empty,
  # e.g. 'make DEBUG=Y'.
  else
    ifneq (,$(findstring ifort,$(FC)))
      FFLAGS := -O0 -g -traceback -check bounds
    else ifneq (,$(findstring gfortran,$(FC)))
      FFLAGS := -O0 -g -fbacktrace -fbounds-check -fno-range-check
    else
      $(info WARNING: No default compilation flags for compiler $(FC).)
    endif
  endif
else
  FFLAGS = $(FOPT)
endif

# Ensure required compiler flags.
ifneq (,$(findstring gfortran,$(FC)))
  # disable compile-time range checking due to large integer
  # constants in math routines (allows for INTEGER*8)
  ifeq (,$(findstring -fno-range-check,$(FFLAGS)))
    # If FFLAGS passed as an argument, must override
    ifeq ("$(origin FFLAGS)","command line")
      override FFLAGS += -fno-range-check
    else
      FFLAGS += -fno-range-check
    endif
  endif
endif

# Flags used for shared object library.
ifeq (,$(findstring -fPIC,$(FFLAGS)))
  FFLAGS += -fPIC
endif

# Flags used to specify DDCalc folder
FFLAGS += -D DDCALC_DIR=\"$(DDCALC_DIR)\" 

# C++ compiler and flags.  Only used for test program.
# CXX_FLIBS are libraries necessary for linking fortran
# routines (compiler specific).
ifneq (,$(findstring "$(origin CXX)","undefined" "default"))
  ifneq (,$(findstring ifort,$(FC)))
    CXX := icc
  else ifneq (,$(findstring gfortran,$(FC)))
    CXX := g++
  endif
endif
ifneq (,$(findstring "$(origin CXXFLAGS)","undefined" "default"))
  ifneq (,$(findstring icc,$(CXX)))
    CXXFLAGS := 
  else ifneq (,$(findstring g++,$(CXX)))
    CXXFLAGS := 
  endif
endif
ifneq (,$(findstring icc,$(CXX)))
  CXX_FLIBS := -lifcore
else ifneq (,$(findstring g++,$(CXX)))
  CXX_FLIBS := -lgfortran 
else
  $(info WARNING: Unsupported C++ compiler $(CXX).  Build may fail due to)
  $(info unknown library dependencies.)
  CXX_FLIBS := 
endif

# Debugging
#$(info FC      =$(FC) ($(origin FC)))
#$(info FFLAGS  =$(FFLAGS) ($(origin FFLAGS)))
#$(info FOPT  =$(FOPT) ($(origin FOPT)))


#################### FILES / TARGETS #######################

# Main programs
fprograms := DDCalc_run DDTest

# Example/test programs
ftestprograms := DDCalc_exampleF
ctestprograms := DDCalc_exampleC

# Fortran sources
fsrc := DDConstants.f90 DDTypes.f90 DDCommandLine.f90 DDNuclear.f90 \
        DDUtils.f90 DDNumerical.f90 DDStats.f90 DDCouplings.f90 \
        DDWIMP.f90 DDInput.f90 DDHalo.f90 DDRates.f90 DDOutput.f90 \
        DDDetectors.f90 DDTabulation.f90                                      
fsrc_final := DDExperiments.f90 DDExHelp.f90 DDCalc.f90
analyses := LUX_2016.f90 LUX_2013.f90 LUX_2015.f90 DARWIN_Ar.f90 \
            DARWIN_Xe.f90 XENON100_2012.f90 PandaX_2016.f90 \
            PICO_2L.f90 PICO_60_F.f90 PICO_60_I.f90 \
            SIMPLE_2014.f90 SuperCDMS_2014.f90

# Include files
fincludes := 
fincludes := $(patsubst %,$(INCLUDE)/%,$(fincludes))
cincludes := DDCalc.hpp DDExperiments.hpp
cincludes := $(patsubst %,$(INCLUDE)/%,$(cincludes)) \

# Libraries
libname   := DDCalc
statlib   := lib$(libname).a
sharedlib := lib$(libname).so
libraries := $(statlib) $(sharedlib)

# Source files.
fsources := $(fincludes) \
            $(patsubst %,$(SRC)/%,$(fsrc)) \
            $(patsubst %,$(SRC)/analyses/%,$(analyses)) \
            $(patsubst %,$(SRC)/%,$(fsrc_final)) \
            $(patsubst %,$(SRC)/%.f90,$(fprograms)) \
            $(patsubst %,$(EXAMPLES)/%.f90,$(ftestprograms))
csources := $(cincludes) \
            $(patsubst %,$(EXAMPLES)/%.cpp,$(ctestprograms))
sources  := $(fsources) $(csources)
fobjects := $(patsubst %.f90,$(BUILD)/%.o,$(fsrc) $(analyses) $(fsrc_final))

# Additional files to include in distribution tar file
extrafiles := Makefile README LICENSE DDCalc.bib

# Distribution tar file and subdirectory
distfile := DDCalc-$(DDDCALC_VERSION).tar.gz
DIST_DIR := DDCalc-$(DDDCALC_VERSION)

# Files to remove when cleaning.
# Note some compilers leave .dSYM files/directories for debugging.
cleanfiles := $(BUILD)/* \
              $(fprograms) $(fprograms:=.dSYM) \
							$(ftestprograms) $(ftestprograms:=.dSYM) \
							$(ctestprograms)                
# ...also tar file and libs
distcleanfiles := $(cleanfiles) $(LIB)/* $(distfile)

#################### DEPENDENCIES ##########################

# General include file dependencies
$(ctestprograms): $(cincludes)

# Some objects and programs have additional include files
#file.o: file_com.h

#################### RULES #################################

# Default: build only the programs
.DEFAULT_GOAL := bin

# Phony targets
.PHONY: all bin lib examples FORCE clean distclean 

# Build everything (except archives)
all: lib bin examples

# Binaries (programs)
bin: $(fprograms)

# Libraries
lib: $(libraries)

# Example/test binaries (programs)
examples: $(ftestprograms) $(ctestprograms)

# Define a do-nothing rule for dummy targets
# (avoids "nothing to be done" messages)
ifort gfortran:
	@:

# Rule for building programs
$(fprograms): $(statlib)
$(fprograms): % : $(SRC)/%.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -I$(BUILD) -o $@ $< $(LIB)/$(statlib) 

# Rule for building test programs (F)
$(ftestprograms): $(statlib)
$(ftestprograms): % : $(EXAMPLES)/%.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -I$(BUILD) -o $@ $< $(LIB)/$(statlib) 

# Rule for building test programs (C++)
$(ctestprograms): $(statlib)
$(ctestprograms): % : $(EXAMPLES)/%.cpp
	@echo "===== Compiling $@ ====="
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -o $@ $< $(LIB)/$(statlib) $(CXX_FLIBS)

# Rules for building objects
$(BUILD)/%.o : $(SRC)/%.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -o $@ -c $<
$(BUILD)/%.o : $(SRC)/analyses/%.f90
	@echo "===== Compiling $@ ====="
	$(FC) $(FFLAGS) -o $@ -c $<

# Rule for building static library
$(statlib): $(fobjects)
	@echo "===== Generating $@ ====="
	ar rcs $(LIB)/$@ $^

# Rule for building shared library
$(sharedlib): $(fobjects)
	@echo "===== Generating $@ ====="
	$(FC) $(FFLAGS) -shared -o $(LIB)/$@ $^

# Generate archive file
dist: $(distfile)
$(distfile): $(sources) $(extrafiles)
	@echo "===== Generating $@ ====="
	mkdir -p ".archive-temp/$(DIST_DIR)/build"
	mkdir -p ".archive-temp/$(DIST_DIR)/lib"
	cp -r -al $(SRC) $(INCLUDE) $(DOC) $(EXAMPLES) \
	 $(DATA) $(extrafiles) ".archive-temp/$(DIST_DIR)/"
	tar -czf $@ -C ".archive-temp" "$(DIST_DIR)"
	rm -rf ".archive-temp"

# Remove files
clean: FORCE
	@echo "===== Cleaning ====="
	-rm -rf $(cleanfiles)

# Remove more files
distclean:FORCE
	@echo "===== Cleaning ====="
	-rm -rf $(distcleanfiles)
