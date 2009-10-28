# Kelvin
# Copyright 2009, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Compiled executables and scripts will be installed in $BINDIR
BINDIR=/usr/local/bin
OWNER=root
GROUP=root

# $INCDIR and $LIBDIR should point to where headers and libraries for
# GSL (GNU Scientific Library) can be found. Remember you can specify these
# as command-line macros, e.g. at OSC:
# $ make INCDIR=/home/ccri0005/include LIBDIR=/home/ccri0005/lib
#
ifndef LIBDIR
LIBDIR=/usr/local/lib
endif
ifndef INCDIR
INCDIR=/usr/local/include
endif
KVNLIBDIR := $(shell pwd)/lib
KELVIN_ROOT := $(shell pwd)
TEST_KELVIN := $(KELVIN_ROOT)/kelvin
TEST_UPDATE := $(KELVIN_ROOT)/seq_update/calc_updated_ppl
VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
PLATFORM_NAME := $(shell echo `uname -m`-`uname -s`)
empty:=
space:= $(empty) $(empty)
PLATFORM = $(subst $(space),-,$(PLATFORM_NAME))
INCFLAGS := -I$(INCDIR)

CC := gcc
#CC := icc # For the Intel C Compiler at OSC
GCCOPT := 2 # GCC optimization level, 0=none, 1=default, 2=some (OSC's recommendation), 3=all
CFLAGS := -Wall -Werror -DGCCOPT=$(GCCOPT) -O$(GCCOPT) # -Wshadow # PitA gcc won't tell me optimization level
CFLAGS += -D_REENTRANT # Thead-safe (different prototype) version of strtok_r under Solaris when using pthread
LDFLAGS := -L$(LIBDIR) -L$(KVNLIBDIR)

# For further details on compilation-time conditionals, see kelvin.c or the Doxygen documentation.

#CFLAGS += -g # Only an ~10% drag on performance and we can monitor running processes w/symbols.
CFLAGS += -fopenmp # Uncomment for multi-threading if using GCC 4.2+. MUST USE GSL TOO.
#CFLAGS += -openmp # Same as above, but only for Intel C Compiler
#LPTM3FLAG = -lptmalloc3 # For ptmalloc3 allocator, some performance gains, tighter memory use w/OpenMP, but not on Mac.
CFLAGS += -DSIMPLEPROGRESS # Simplify progress reporting to a wobbly percentage and estimated time left
#CFLAGS += -DMEMSTATUS # Display time and memory consumption every 30 seconds
#CFLAGS += -DMEMGRAPH # Log terse time and memory consumption info to a data file every 30 seconds for graphing
#CFLAGS += -DPOLYSTATISTICS # Display extensive polynomial statistics every 2Mp and at milestones
#CFLAGS += -DDMUSE # For our own static memory management, not beneficial as yet.
#CFLAGS += -DDMTRACK # For our own memory tracking
#CFLAGS += -DTREEEVALUATE # Use evaluateValue of tree instead of evaluatePoly of list.
#CFLAGS += -DFAKEEVALUATE # Don't evaluate at all - use only for exercise/compilation. Results will be wrong!
#CFLAGS += -DPOLYUSE_DL # Dynamically load compiled polynomials for in-process use
#ADD_LDFLAGS += -ldl # ditto
#CFLAGS += -DPOLYCODE_DL # Enable generation of dynamic library code for selected polynomials
#CFLAGS += -DPOLYCOMP_DL # Enable compilation of dynamic library code for selected polynomials
#CFLAGS += -DPOLYCHECK_DL # Keep both built and compiled DL polys and compare results (can be noisy!)
#CFLAGS += -DTELLRITA # Relay all log messages to rita via UDP
#ADD_LDFLAGS += -lsocket -lnsl # ditto for under Solaris
#CFLAGS += -DUSE_SSD # Experimental use of solid state drive when building polynomials. NOT THREAD-SAFE!
CFLAGS += -DUSE_GSL # Use GNU Scientific Library (GSL) statistical routines instead of internal ones
ADD_LDFLAGS += -lgsl -lgslcblas -lm # ditto
#CFLAGS += -DVERIFY_GSL # Use both internal and GSL returning internal and printing if error > 1e-13, no OpenMP

LDFLAGS += ${ADD_LDFLAGS}
export KVNLIBDIR VERSION CC CFLAGS LDFLAGS INCFLAGS KELVIN_ROOT TEST_KELVIN

KOBJS = kelvin.o dcuhre.o
OBJS = kelvinInit.o kelvinTerm.o iterationSupport.o integrationSupport.o \
	kelvinHandlers.o kelvinWriteFiles.o dkelvinWriteFiles.o \
	ppl.o saveResults.o trackProgress.o \
	summary_result.o tp_result_hash.o

INCS = kelvin.h kelvinGlobals.h kelvinLocals.h kelvinHandlers.h \
	kelvinInit.h kelvinTerm.h \
	iterationLocals.h iterationSupport.h \
	integrationGlobals.h integrationLocals.h integrationSupport.h \
	kelvinWriteFiles.h dkelvinWriteFiles.h \
	ppl.h dcuhre.h saveResults.h summary_result.h trackProgress.h tp_result_hash.h

# Binary releases include kelvin_$(PLATFORM)
all : kelvin seq_update/calc_updated_ppl 

install : $(BINDIR)/kelvin-$(VERSION) \
          $(BINDIR)/calc_updated_ppl \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh

kelvin : libs $(KOBJS) $(OBJS) $(INCS)
	$(CC) -o $@ $(KOBJS) $(OBJS) -lped -lconfig -lklvnutls -lm -lpthread $(LDFLAGS) $(CFLAGS) $(EXTRAFLAG)

kelvin_$(PLATFORM) : libs $(KOBJS) $(OBJS)
	$(CC) -static $(LPTMFLAG) -o $@ $(KOBJS) $(OBJS) $(LDFLAGS) $(CFLAGS) $(EXTRAFLAG)

.PHONY : seq_update/calc_updated_ppl
seq_update/calc_updated_ppl :
	+make -C seq_update -f Makefile calc_updated_ppl

%.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) $< -o $@

.PHONY : libs
libs :
	+make -C utils -f Makefile all
	+make -C config -f Makefile all
	+make -C pedlib -f Makefile all

.PHONY : clean
clean :
	make -C pedlib -f Makefile clean
	make -C config -f Makefile clean
	make -C utils -f Makefile clean
	make -C seq_update -f Makefile clean
	rm -f $(KOBJS) $(OBJS) kelvin seq_update/calc_updated_ppl lib/*
	make -C test-suite -f Makefile clean

.PHONY : test test-USE_DL
test-USE_DL :
	+make -C test-suite -f Makefile test-USE_DL

.PHONY : test test-FIXED
test-FIXED :
	+make -C test-suite -f Makefile test-FIXED
test :
	+make -C test-suite -f Makefile test

$(BINDIR)/kelvin-$(VERSION) : kelvin
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin $(BINDIR)/kelvin-$(VERSION)

$(BINDIR)/calc_updated_ppl : seq_update/calc_updated_ppl
	install -o $(OWNER) -g $(GROUP) -m 0755 -p seq_update/calc_updated_ppl $(BINDIR)/calc_updated_ppl

$(BINDIR)/convert_br.pl : seq_update/convert_br.pl
	install -o $(OWNER) -g $(GROUP) -m 0755 -p seq_update/convert_br.pl $(BINDIR)/convert_br.pl

$(BINDIR)/compileDL.sh : compileDL.sh
	install -o $(OWNER) -g $(GROUP) -m 0755 -p compileDL.sh $(BINDIR)/compileDL.sh
