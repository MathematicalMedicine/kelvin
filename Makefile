# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Compiled executables and scripts will be installed in $BINDIR
BINDIR=/usr/local/bin

# $INCDIR and $LIBDIR should point to where headers and libraries for
# GSL (GNU Scientific Library) can be found
INCDIR=/usr/local/include
LIBDIR=/usr/local/lib
KVNLIBDIR := $(shell pwd)/lib
KVNINCDIR := $(shell pwd)/include
VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
INCFLAGS := -I$(INCDIR) -I$(KVNINCDIR)

CC := gcc
CFLAGS := -Wall -O3
LDFLAGS := -L$(LIBDIR) -L$(KVNLIBDIR) -lped -lutils -lgsl -lgslcblas

CFLAGS += -fopenmp # Uncomment BOTH of these if you have an OpenMP-capable compiler...
LDFLAGS += -fopenmp # ...and want to use multiple threads for evaluations.

export KVNLIBDIR KVNINCDIR VERSION CC CFLAGS LDFLAGS INCFLAGS

KOBJS = kelvin.o
DKOBJS = dkelvin.o dcuhre.o
OBJS = ppl.o config.o saveResults.o
INCS = kelvin.h dkelvin.h dcuhre.h saveResults.h

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

$(BINDIR)/kelvin-$(VERSION) : kelvin
	install -o root -g root -m 0755 -p kelvin $(BINDIR)/kelvin-$(VERSION)

$(BINDIR)/dkelvin-$(VERSION) : dkelvin
	install -o root -g root -m 0755 -p dkelvin $(BINDIR)/dkelvin-$(VERSION)

$(BINDIR)/calc_updated_ppl : calc_updated_ppl
	install -o root -g root -m 0755 -p calc_updated_ppl $(BINDIR)/calc_updated_ppl

$(BINDIR)/seq_update_avghet.pl : seq_update/seq_update_avghet.pl
	install -o root -g root -m 0755 -p seq_update/seq_update_avghet.pl $(BINDIR)/seq_update_avghet.pl



