# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Compiled executables and scripts will be installed in $BINDIR
BINDIR=/usr/local/bin

# $INCDIR and $LIBDIR should point to where headers and libraries for
# GSL (GNU Scientific Library) can be found. Remember you can specify these
# as command-line macros, e.g. at OSC:
# $ make INCDIR=/home/ccri0005/include LIBDIR=/home/ccri0005/lib
#
INCDIR=/usr/local/include
LIBDIR=/usr/local/lib
KVNLIBDIR := $(shell pwd)/lib
KVNINCDIR := $(shell pwd)/include
VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
PLATFORM_NAME := $(shell echo `uname -m`-`uname -s`)
empty:=
space:= $(empty) $(empty)
PLATFORM = $(subst $(space),-,$(PLATFORM_NAME))
INCFLAGS := -I$(INCDIR) -I$(KVNINCDIR)

CC := gcc
#CC := icc # For the Intel C Compiler at OSC
CFLAGS := -Wall # -O3 -Wshadow

# We currently do not use GSL, so use this line:
LDFLAGS := -L$(LIBDIR) -L$(KVNLIBDIR) -lped -lkutils -lm -lpthread
# Alternatively, to use or test with GSL...
#CFLAGS += -DUSE_GSL
#LDFLAGS := -L$(LIBDIR) -L$(KVNLIBDIR) -lped -lkutils -lgsl -lgslcblas -lm -lpthread

# For further details on compilation-time conditionals, see kelvin.c or the Doxygen documentation.

#CFLAGS += -g # Only an ~10% drag on performance and we can monitor running processes w/symbols.
CFLAGS += -fopenmp # Uncomment if you have an OpenMP-capable compiler and want to use multiple threads for evaluations.
#CFLAGS += -openmp # Same as above, but only for Intel C Compiler
LPTM3FLAG = -lptmalloc3 # For ptmalloc3 allocator, some performance gains, tighter memory use w/OpenMP, but not on Mac.
CFLAGS += -DSIMPLEPROGRESS # Simplify progress reporting to a wobbly percentage and estimated time left
#CFLAGS += -DMEMSTATUS # Display time and memory consumption every 30 seconds
#CFLAGS += -DMEMGRAPH # Log terse time and memory consumption info to a data file every 30 seconds for graphing
#CFLAGS += -DPOLYSTATISTICS # Display extensive polynomial statistics every 2Mp and at milestones
#CFLAGS += -DDMUSE # For our own static memory management, not beneficial as yet.
#CFLAGS += -DDMTRACK # For our own memory tracking
#CFLAGS += -DTREEEVALUATE # Use evaluateValue of tree instead of evaluatePoly of list.
#CFLAGS += -DFAKEEVALUATE # Don't evaluate at all - use only for exercise/compilation. Results will be wrong!
#CFLAGS += -DPOLYUSE_DL # Dynamically load compiled polynomials for in-process use
#ADD_LDFLAGS += -ldl # "
#CFLAGS += -DPOLYCODE_DL # Enable generation of dynamic library code for selected polynomials
#CFLAGS += -DPOLYCOMP_DL # Enable compilation of dynamic library code for selected polynomials
#CFLAGS += -DPOLYCHECK_DL # Keep both built and compiled DL polys and compare results (can be noisy!)
#CFLAGS += -DTELLRITA # Relay all log messages to rita via UDP
#CFLAGS += -DUSE_SSD # Experimental use of solid state drive when building polynomials. NOT THREAD-SAFE!

LDFLAGS += ${ADD_LDFLAGS}
export KVNLIBDIR KVNINCDIR VERSION CC CFLAGS LDFLAGS INCFLAGS

KOBJS = kelvin.o dcuhre.o
OBJS = ppl.o config.o saveResults.o trackProgress.o kelvinHandlers.o
INCS = kelvin.h dcuhre.h saveResults.h trackProgress.h kelvinHandlers.h \
	kelvinGlobals.h iterationGlobals.h integrationGlobals.h \
	kelvinLocals.h iterationLocals.h integrationLocals.h \
	kelvinInit.c kelvinTerm.c iterationMain.c integrationMain.c kelvinWriteFiles.c

# Binary releases include kelvin_$(PLATFORM)
all : kelvin calc_updated_ppl 

install : $(BINDIR)/kelvin-$(VERSION) \
          $(BINDIR)/calc_updated_ppl \
          $(BINDIR)/seq_update_br.pl \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh

kelvin : libs $(KOBJS) $(OBJS)
	$(CC) -o $@ $(KOBJS) $(OBJS) $(LDFLAGS) $(CFLAGS) $(EXTRAFLAG) $(LPTMFLAG)

kelvin_$(PLATFORM) : libs $(KOBJS) $(OBJS)
	$(CC) -static $(LPTMFLAG) -o $@ $(KOBJS) $(OBJS) $(LDFLAGS) $(CFLAGS) $(EXTRAFLAG)

calc_updated_ppl : seq_update/calc_updated_ppl.c
	$(CC) -o $@ $(CFLAGS) $(INCFLAGS) seq_update/calc_updated_ppl.c -lm

%.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) $< -o $@

.PHONY : libs
libs :
	make -C utils -f Makefile all
	make -C pedlib -f Makefile all

.PHONY : clean
clean :
	make -C pedlib -f Makefile clean
	make -C utils -f Makefile clean
	rm -f $(KOBJS) $(OBJS) kelvin calc_updated_ppl
	make -C test-suite -f Makefile clean

.PHONY : test test-USE_DL
test-USE_DL :
	make -C test-suite -f Makefile test-USE_DL
test :
	make -C test-suite -f Makefile test

$(BINDIR)/kelvin-$(VERSION) : kelvin
	install -o root -g root -m 0755 -p kelvin $(BINDIR)/kelvin-$(VERSION)

$(BINDIR)/calc_updated_ppl : calc_updated_ppl
	install -o root -g root -m 0755 -p calc_updated_ppl $(BINDIR)/calc_updated_ppl

$(BINDIR)/seq_update_br.pl : seq_update/seq_update_br.pl
	install -o root -g root -m 0755 -p seq_update/seq_update_br.pl $(BINDIR)/seq_update_br.pl

$(BINDIR)/convert_br.pl : seq_update/convert_br.pl
	install -o root -g root -m 0755 -p seq_update/convert_br.pl $(BINDIR)/convert_br.pl

$(BINDIR)/compileDL.sh : compileDL.sh
	install -o root -g root -m 0755 -p compileDL.sh $(BINDIR)/compileDL.sh
