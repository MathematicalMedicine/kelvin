# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Compiled executables and scripts will be installed in $BINDIR
BINDIR=/usr/local/bin

# $INCDIR and $LIBDIR should point to where headers and libraries for
# GSL (GNU Scientific Library) can be found. Remember you can specify these
# as command-line macros, e.g. $ make "INCDIR=$HOME/usr/local/include"
INCDIR=/usr/local/include
LIBDIR=/usr/local/lib
KVNLIBDIR := $(shell pwd)/lib
KVNINCDIR := $(shell pwd)/include
VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
INCFLAGS := -I$(INCDIR) -I$(KVNINCDIR)

CC := gcc
CFLAGS := -Wall -O3
LDFLAGS := -L$(LIBDIR) -L$(KVNLIBDIR) -lped -lutils -lgsl -lgslcblas

# For further details on compilation-time conditionals, see kelvin.c or the Doxygen documentation.

#CFLAGS += -g # Only an ~10% drag on performance and we can monitor running processes w/symbols.
CFLAGS += -fopenmp # Uncomment if you have an OpenMP-capable compiler and want to use multiple threads for evaluations.
#LDFLAGS += -lptmalloc3 # For ptmalloc3 allocator, some performance gains, tighter memory use w/OpenMP, but not on Mac.
#CFLAGS += -DSIMPLEPROGRESS # Simplify progress reporting to a wobbly percentage and estimated time left
#CFLAGS += -DMEMSTATUS # Display time and memory consumption every 30 seconds
#CFLAGS += -DMEMGRAPH # Log terse time and memory consumption info to a data file every 30 seconds for graphing
#CFLAGS += -DPOLYSTATISTICS # Display extensive polynomial statistics every 2Mp and at milestones
#CFLAGS += -DDMUSE # For our own static memory management, not beneficial as yet.
#CFLAGS += -DDMTRACK # For our own memory tracking
#CFLAGS += -DTREEEVALUATE # Use evaluateValue of tree instead of evaluatePoly of list.
#CFLAGS += -DFAKEEVALUATE # Don't evaluate at all - use only for exercising build. Results will be wrong!
#CFLAGS += -DPOLYCOMP # Enable compilation and distribution of selected polynomials for faster evaluation
#CFLAGS += -DPOLYCOMP_DL -ldl # Dynamically load compiled polynomials for in-process use

export KVNLIBDIR KVNINCDIR VERSION CC CFLAGS LDFLAGS INCFLAGS

KOBJS = kelvin.o
DKOBJS = dkelvin.o dcuhre.o
OBJS = ppl.o config.o saveResults.o trackProgress.o
INCS = kelvin.h dkelvin.h dcuhre.h saveResults.h trackProgress.h

all : kelvin dkelvin calc_updated_ppl

install : $(BINDIR)/kelvin-$(VERSION) \
          $(BINDIR)/dkelvin-$(VERSION) \
          $(BINDIR)/calc_updated_ppl \
          $(BINDIR)/seq_update_avghet.pl

kelvin : libs $(KOBJS) $(OBJS)
	$(CC) -o $@ $(KOBJS) $(OBJS) $(LDFLAGS) $(CFLAGS)

dkelvin : libs $(DKOBJS) $(OBJS)
	$(CC) -o $@ $(DKOBJS) $(OBJS) $(LDFLAGS) $(CFLAGS)

calc_updated_ppl : seq_update/calc_updated_ppl.c
	$(CC) -o $@ $(CFLAGS) seq_update/calc_updated_ppl.c

%.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : libs
libs :
	make -C utils -f Makefile all
	make -C pedlib -f Makefile all

.PHONY : clean
clean :
	make -C pedlib -f Makefile clean
	make -C utils -f Makefile clean
	rm -f $(KOBJS) $(DKOBJS) $(OBJS) kelvin dkelvin calc_updated_ppl

.PHONY : test
test :
	make -C test-suite -f Makefile all

$(BINDIR)/kelvin-$(VERSION) : kelvin
	install -o root -g root -m 0755 -p kelvin $(BINDIR)/kelvin-$(VERSION)

$(BINDIR)/dkelvin-$(VERSION) : dkelvin
	install -o root -g root -m 0755 -p dkelvin $(BINDIR)/dkelvin-$(VERSION)

$(BINDIR)/calc_updated_ppl : calc_updated_ppl
	install -o root -g root -m 0755 -p calc_updated_ppl $(BINDIR)/calc_updated_ppl

$(BINDIR)/seq_update_avghet.pl : seq_update/seq_update_avghet.pl
	install -o root -g root -m 0755 -p seq_update/seq_update_avghet.pl $(BINDIR)/seq_update_avghet.pl



