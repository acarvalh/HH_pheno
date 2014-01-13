################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all

GCC = g++

# main programs and required modules 

#MAIN =  4jet-reco-signal  #4jet-reco-signal-test #4jet-reco-back  4jet-reco-signal-gamma 4jet-reco-back-gamma
MAIN =  HH_VBF  #4jet-reco-signal-test #4jet-reco-back  4jet-reco-signal-gamma 4jet-reco-back-gamma
MODULES = Functions #

# search path for modules

MDIR = ./

# root
ROOTINCS = $(shell /home/xanda/root/bin/root-config --cflags) 
ROOTLIBS = $(shell /home/xanda/root/bin/root-config --glibs) 

#LHAPDF
#LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
#LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
#LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# fastjet
FJINCS = $(shell /home/xanda/Documents/fastjet-install/bin/fastjet-config --cxxflags)
FJCLIBS    = $(shell /home/xanda/Documents/fastjet-install/bin/fastjet-config --libs)

# additional include directories
INCPATH = -I../include $(LHAPDFINCS) $(FJINCS) $(ROOTINCS)

# additional libraries to be included 
LIBS = $(LHAPDFLIBS) $(FJCLIBS)   $(ROOTLIBS)

# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -Wall -ansi -O3 

############################## do not change ###################################

SHELL=/bin/bash

CC=$(GCC)

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(INCPATH)

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(LIBS)

-include $(addsuffix .d,$(PGMS))

# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.cc Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.cc Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
	echo "delete old executables"           


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory 

clean:
	@ -rm -rf *.d *.o *~ $(MAIN) *.eps *.data plots/*~ plots/*.eps *.a analysis/*~ analysis/*.eps
.PHONY: clean

lib:
	ar rcs libResPairTagger.a  ResPairTagger.o 

################################################################################
