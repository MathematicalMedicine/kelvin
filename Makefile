# Kelvin
# Copyright 2009, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes
# only.

## Variables and options in this Makefile are fully documented in
## doc/compileoptions.html

## Directory into which compiled executables and scripts will be installed.
BINDIR=/usr/local/bin

## User and group IDs by which installed execuatbles and scripts will be owned.
OWNER=root
GROUP=root

## Directories in which optional header files and libraries can be found (GSL,
## etc). Remember you can specify these as command-line macros, e.g. at OSC:
## $ make INCDIR=/home/ccri0005/include LIBDIR=/home/ccri0005/lib
ifndef LIBDIR
  LIBDIR=/usr/local/lib
endif
ifndef INCDIR
  INCDIR=/usr/local/include
endif

## The C compiler to be used to build executables. Pick one.
## GCC (GNU C Compiler)
CC := gcc
## ICC (Intel C Compiler)
# CC := icc

## GCC optimization level, 0=none, 1=default, 2=some (recommended), 3=all
GCCOPT := 2

## Enable OpenMP support. Requires icc or gcc 4.2+, and GSL
USE_OPENMP := no

## Enable use of GSL (GNU Scientific Library). Don't forget to set
## INCDIR and LIBDIR (above) accordingly.
USE_GSL := yes

## Enable use of ptmalloc3. Don't forget to set LIBDIR (above) accordingly.
## Not available on OSX.
# USE_PTMALLOC3 := yes

## Enable use of Hoard. Don't forget to set LIBDIR (above) accordingly.
# USE_HOARD := yes

## Beginning of Kelvin-specific options
CFLAGS :=

##                                                     ##
## Should be no need to make changes beyond this point ##
##                                                     ##

LDFLAGS := -dynamic 
ADD_LDFLAGS :=

CFLAGS += -O$(GCCOPT)
CFLAGS += -DGCCOPT=$(GCCOPT)

## Compiler warnings. In development, abort on warning
CFLAGS += -Wall
# CFLAGS += -Wshadow -Werror

## Enable debugging symbols. Inflicts a small drag (10%) on performance
CFLAGS += -g

# Required to get re-entrant routines under Solaris, benign on other platforms
CFLAGS += -D_REENTRANT


# If OpenMP support has been enabled, GSL is required. The GSL-replacement 
# routines are not thread-safe.
ifeq ($(strip $(USE_OPENMP)), yes)
  USE_GSL := yes
  ifeq ($(strip $(CC)), gcc)
    # Compiler flags for GCC
    CFLAGS += -fopenmp
    USE_PTHREAD := yes
    ADD_LDFLAGS += -fopenmp
  else
    ifeq ($(strip $(CC)), icc)
      # Compiler flags for ICC
      CFLAGS += -openmp 
      USE_PTHREAD := yes
      ADD_LDFLAGS += -openmp
    endif
  endif
endif

# If GSL support has been enabled
ifeq ($(strip $(USE_GSL)), yes)
  CFLAGS += -DUSE_GSL
  ADD_LDFLAGS += -lgsl -lgslcblas -lm
endif

# If ptmalloc3 support has been enabled
ifeq ($(strip $(USE_PTMALLOC3)), yes)
  ADD_LDFLAGS += -lptmalloc3
  USE_PTHREAD := yes
endif

# If Hoard support has been enabled
ifeq ($(strip $(USE_HOARD)), yes)
  ADD_LDFLAGS += -lhoard
endif

# If libpthread is required, either for OpenMP or libptmalloc3
ifeq ($(strip $(USE_PTHREAD)), yes)
  ADD_LDFLAGS += -lpthread
endif

PLATFORM_NAME := $(shell echo `uname -m`-`uname -s`)
empty:=
space:= $(empty) $(empty)
PLATFORM = $(subst $(space),$(empty),$(PLATFORM_NAME))
KVNLIBDIR := $(shell pwd)/lib
KELVIN_ROOT := $(shell pwd)
TEST_KELVIN := $(KELVIN_ROOT)/kelvin.$(PLATFORM)
TEST_UPDATE := $(KELVIN_ROOT)/seq_update/calc_updated_ppl

# If we're building in an svn-managed context, get AND preserve the latest svn version
SVNVERSION := $(subst exported,,$(shell svnversion 2>/dev/null))
ifeq ("$(strip $(SVNVERSION))","")
  SVNVERSION := $(shell cat .svnversion)
else
  UPDATE_SVNVERSION := $(shell echo $(SVNVERSION) >.svnversion)
endif
VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)

INCFLAGS := -I$(INCDIR)

# testmac doesn't recognize the -rdynamic bit...
LDFLAGS := -rdynamic -L$(LIBDIR) -L$(KVNLIBDIR)

# Flags for BCMM use only

#CFLAGS += -DDISTRIBUTION # Eliminates all diagnostics for distribution purposes
ifneq (,$(wildcard /usr/include/execinfo.h))
#  CFLAGS += -DBACKTRACE # Add backtrace where supported
endif
#CFLAGS += -DMEMSTATUS # Display time and memory consumption every 30 seconds
#CFLAGS += -DMEMGRAPH # Log terse time and memory consumption info to a data file every 30 seconds for graphing
#CFLAGS += -DPOLYSTATISTICS # Display extensive polynomial statistics every raw 8Mp and at milestones
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
#CFLAGS += -DFULLLOG # Write all log messages to kelvin.full_log if TELLRITA isn't working
#ADD_LDFLAGS += -lsocket -lnsl # ditto for under Solaris
#CFLAGS += -DUSE_SSD # Experimental use of solid state drive when building polynomials. NOT THREAD-SAFE!
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

# Binary releases include kelvin.$(PLATFORM)
all : kelvin-$(VERSION) seq_update/calc_updated_ppl

dist :
	- rm -rf kelvin-$(SVNVERSION)
	mkdir kelvin-$(SVNVERSION)
	mkdir kelvin-$(SVNVERSION)/bin
	ln bin/kelvin.* kelvin-$(SVNVERSION)/bin
	ln .maj .min .pat .svnversion Kelvin CHANGES COPYRIGHT Makefile PedCount.pl kf.pm convertconfig.pl rebuild.sh *.[ch] kelvin-$(SVNVERSION)
	mkdir kelvin-$(SVNVERSION)/lib
	mkdir kelvin-$(SVNVERSION)/utils
	ln utils/Makefile utils/*.[ch] utils/wordDiff.pl kelvin-$(SVNVERSION)/utils
	mkdir kelvin-$(SVNVERSION)/pedlib
	ln pedlib/Makefile pedlib/*.[ch] kelvin-$(SVNVERSION)/pedlib
	mkdir kelvin-$(SVNVERSION)/config
	ln config/Makefile config/*.[ch] kelvin-$(SVNVERSION)/config
	mkdir kelvin-$(SVNVERSION)/seq_update
	ln seq_update/Makefile seq_update/*.[ch] kelvin-$(SVNVERSION)/seq_update
	mkdir -p kelvin-$(SVNVERSION)/test-suite/dynamic-grid/PE/SA_DT
	ln test-suite/Makefile kelvin-$(SVNVERSION)/test-suite
	ln test-suite/dynamic-grid/Makefile kelvin-$(SVNVERSION)/test-suite/dynamic-grid
	ln test-suite/dynamic-grid/PE/Makefile kelvin-$(SVNVERSION)/test-suite/dynamic-grid/PE
	ln test-suite/dynamic-grid/PE/SA_DT/* kelvin-$(SVNVERSION)/test-suite/dynamic-grid/PE/SA_DT
	mkdir -p kelvin-$(SVNVERSION)/test-suite/PedCount/No_config
	ln test-suite/PedCount/Makefile kelvin-$(SVNVERSION)/test-suite/PedCount
	ln test-suite/PedCount/No_config/* kelvin-$(SVNVERSION)/test-suite/PedCount/No_config
	mkdir -p kelvin-$(SVNVERSION)/test-suite/seq_update/d-2pt-le
	ln test-suite/seq_update/Makefile kelvin-$(SVNVERSION)/test-suite/seq_update
	ln test-suite/seq_update/d-2pt-le/* kelvin-$(SVNVERSION)/test-suite/seq_update/d-2pt-le
	mkdir kelvin-$(SVNVERSION)/doc
	ln doc/*.html doc/*.png doc/*.gif kelvin-$(SVNVERSION)/doc
	tar -hcvzf kelvin-$(SVNVERSION).tar.gz kelvin-$(SVNVERSION)/
#	rm -rf kelvin-$(SVNVERSION)

install : $(BINDIR)/kelvin-$(VERSION) \
          $(BINDIR)/calc_updated_ppl \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh \
	  $(BINDIR)/PedCount.pl \
	  $(BINDIR)/kf.pm \
	  $(BINDIR)/Kelvin

install-prebuilt : $(BINDIR)/calc_updated_ppl \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh \
	  $(BINDIR)/PedCount.pl \
	  $(BINDIR)/kf.pm \
	  $(BINDIR)/Kelvin
	install -o $(OWNER) -g $(GROUP) -m 0755 -p bin/kelvin.$(PLATFORM) $(BINDIR)/kelvin-$(VERSION)

.PHONY : kelvin
kelvin : kelvin-$(VERSION)

kelvin-$(VERSION) : libs $(KOBJS) $(OBJS) $(INCS)
	$(CC) -o $@ $(KOBJS) $(OBJS) -lped -lconfig -lklvnutls -lm -lpthread $(LDFLAGS) $(CFLAGS) $(EXTRAFLAG)
	cp $@ $@-$(SVNVERSION)
	cp $@ bin/kelvin.$(PLATFORM)

.PHONY : seq_update/calc_updated_ppl
seq_update/calc_updated_ppl :
	+make -C seq_update -f Makefile calc_updated_ppl

%.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) -DVERSION='"V$(VERSION)"' -DSVNVERSION='"$(SVNVERSION)"' $< -o $@

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
	rm -f $(KOBJS) $(OBJS) kelvin.$(PLATFORM) seq_update/calc_updated_ppl lib/libconfig.a lib/klvnutls.a lib/libped.a
	make -C test-suite -f Makefile clean

.PHONY : test-USE_DL
# Polynomial compilation tests
test-USE_DL :
	+make -C test-suite -f Makefile test-USE_DL

.PHONY : test-FIXED
# Internal dynamic-to-fixed grid validation tests
test-FIXED :
	+make -C test-suite -f Makefile test-FIXED

.PHONY : test
# Extensive capability tests (~30 minutes)
test :
	+make -C test-suite -f Makefile test

.PHONY : check
# Brief installation validation (~1-2 minutes)
check :
	+make -C test-suite -f Makefile check

$(BINDIR)/kelvin-$(VERSION) : kelvin-$(VERSION)
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION) $(BINDIR)/kelvin-$(VERSION)

$(BINDIR)/calc_updated_ppl : seq_update/calc_updated_ppl
	install -o $(OWNER) -g $(GROUP) -m 0755 -p seq_update/calc_updated_ppl $(BINDIR)/calc_updated_ppl

$(BINDIR)/convert_br.pl : seq_update/convert_br.pl
	install -o $(OWNER) -g $(GROUP) -m 0755 -p seq_update/convert_br.pl $(BINDIR)/convert_br.pl

$(BINDIR)/compileDL.sh : compileDL.sh
	install -o $(OWNER) -g $(GROUP) -m 0755 -p compileDL.sh $(BINDIR)/compileDL.sh

$(BINDIR)/PedCount.pl : PedCount.pl
	install -o $(OWNER) -g $(GROUP) -m 0755 -p PedCount.pl $(BINDIR)/PedCount.pl

$(BINDIR)/kf.pm : kf.pm
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kf.pm $(BINDIR)/kf.pm

$(BINDIR)/Kelvin : Kelvin
	perl -pe "s|NO_KELVIN_ROOT|$(BINDIR)|; s|NO_KELVIN_BINARY|$(BINDIR)/kelvin-$(VERSION)|; s|NO_SEQUPDATE_BINARY|$(BINDIR)/calc_updated_ppl|; s|NO_PEDCOUNT_SCRIPT|$(BINDIR)/PedCount.pl|" Kelvin > Kelvin.local
	install -o $(OWNER) -g $(GROUP) -m 0755 -p Kelvin.local $(BINDIR)/Kelvin
