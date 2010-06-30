# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes
# only.

## Variables and options in this Makefile are fully documented in
## doc/compileoptions.html

## Directory into which compiled executables and scripts will be installed.
BINDIR=~/mykelvin

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
USE_OPENMP := yes

## Enable use of GSL (GNU Scientific Library). Don't forget to set
## INCDIR and LIBDIR (above) accordingly.
USE_GSL := yes

## Enable use of ptmalloc3. Don't forget to set LIBDIR (above) accordingly.
## Not available on OSX.
# USE_PTMALLOC3 := yes

## Enable use of Hoard. Don't forget to set LIBDIR (above) accordingly.
# USE_HOARD := yes

## Enable use of MySQL database for HLOD storage/retrieval/distributed processing.
# USE_STUDYDB := yes

## Beginning of Kelvin-specific options
CFLAGS :=

##                                                     ##
## Should be no need to make changes beyond this point ##
##                                                     ##

ABSBINDIR=$(shell echo $(BINDIR))
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

VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
PLATFORM_NAME := $(shell echo `uname -m`-`uname -s`)
empty:=
space:= $(empty) $(empty)
PLATFORM = $(subst $(space),$(empty),$(PLATFORM_NAME))
KVNLIBDIR := $(shell pwd)/lib
KELVIN_ROOT := $(shell pwd)
TEST_KELVIN := $(KELVIN_ROOT)/kelvin-$(VERSION)
KELVIN_SCRIPT := $(KELVIN_ROOT)/Kelvin
PEDCOUNT_SCRIPT := ($KELVIN_ROOT)/PedCount.pl
SEQUPDATE_BINARY := $(KELVIN_ROOT)/seq_update/calc_updated_ppl

# If we're building in an svn-managed context, get AND preserve the latest svn version
SVNVERSION := $(subst exported,,$(shell svnversion 2>/dev/null))
ifeq ("$(strip $(SVNVERSION))","")
  SVNVERSION := $(shell cat .svnversion)
else
  UPDATE_SVNVERSION := $(shell echo $(SVNVERSION) >.svnversion)
endif

INCFLAGS := -I$(INCDIR)

# testmac doesn't recognize the -rdynamic bit...
LDFLAGS := -rdynamic -L$(LIBDIR) -L$(KVNLIBDIR)

# Flags for BCMM use only

#CFLAGS += -DDISTRIBUTION # Eliminates all diagnostics for distribution, don't change, its a dist search target
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

ifeq ($(strip $(USE_STUDYDB)), yes)
  CFLAGS += -DSTUDYDB -I/usr/local/mysql/include -I/usr/include/mysql
  LDFLAGS += -lklvndb -lmysqlclient -L/usr/local/mysql/lib -L/usr/lib64/mysql
endif

LDFLAGS += ${ADD_LDFLAGS}
export KVNLIBDIR VERSION CC CFLAGS LDFLAGS INCFLAGS KELVIN_ROOT TEST_KELVIN KELVIN_SCRIPT PEDCOUNT_SCRIPT SEQUPDATE_BINARY


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

all : kelvin-$(VERSION) seq_update/calc_updated_ppl

dist :
	- rm -rf kelvin-$(VERSION)
	mkdir kelvin-$(VERSION)
	mkdir kelvin-$(VERSION)/bin
	ln bin/kelvin.* kelvin-$(VERSION)/bin
	ln bin/calc_updated_ppl.* kelvin-$(VERSION)/bin
	ln README .maj .min .pat .svnversion Kelvin CHANGES COPYRIGHT PedCount.pl kf.pm convertconfig.pl rebuild.sh *.[ch] compileDL.sh kelvin-$(VERSION)
	perl -pe "s|#CFLAGS \+\= \-DDISTRIBUTION|CFLAGS \+\= \-DDISTRIBUTION|;" Makefile > kelvin-$(VERSION)/Makefile
	mkdir kelvin-$(VERSION)/lib
	mkdir kelvin-$(VERSION)/utils
	ln utils/Makefile utils/*.[ch] utils/wordDiff.pl kelvin-$(VERSION)/utils
	mkdir kelvin-$(VERSION)/pedlib
	ln pedlib/Makefile pedlib/*.[ch] kelvin-$(VERSION)/pedlib
	mkdir kelvin-$(VERSION)/config
	ln config/Makefile config/*.[ch] kelvin-$(VERSION)/config
	mkdir kelvin-$(VERSION)/seq_update
	ln seq_update/Makefile seq_update/*.[ch] seq_update/*.pl kelvin-$(VERSION)/seq_update
	mkdir -p kelvin-$(VERSION)/test-suite/dynamic-grid/PE/SA_DT
	ln test-suite/Makefile kelvin-$(VERSION)/test-suite
	ln test-suite/dynamic-grid/Makefile kelvin-$(VERSION)/test-suite/dynamic-grid
	ln test-suite/dynamic-grid/PE/Makefile kelvin-$(VERSION)/test-suite/dynamic-grid/PE
	ln test-suite/dynamic-grid/PE/SA_DT/* kelvin-$(VERSION)/test-suite/dynamic-grid/PE/SA_DT
	mkdir -p kelvin-$(VERSION)/test-suite/PedCount/No_config
	ln test-suite/PedCount/Makefile kelvin-$(VERSION)/test-suite/PedCount
	ln test-suite/PedCount/No_config/* kelvin-$(VERSION)/test-suite/PedCount/No_config
	mkdir -p kelvin-$(VERSION)/test-suite/seq_update/d-2pt-le
	ln test-suite/seq_update/Makefile kelvin-$(VERSION)/test-suite/seq_update
	ln test-suite/seq_update/d-2pt-le/* kelvin-$(VERSION)/test-suite/seq_update/d-2pt-le
	mkdir kelvin-$(VERSION)/doc
	ln doc/*.html doc/*.png doc/*.gif doc/*.txt kelvin-$(VERSION)/doc
	tar -hcvzf kelvin-$(VERSION).tar.gz kelvin-$(VERSION)/
	rm -rf kelvin-$(VERSION)

install : $(BINDIR)/kelvin-$(VERSION) \
          $(BINDIR)/calc_updated_ppl \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh \
	  $(BINDIR)/PedCount.pl \
	  $(BINDIR)/kf.pm \
	  $(BINDIR)/Kelvin

install-prebuilt : $(BINDIR)/kelvin.$(PLATFORM) \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh \
	  $(BINDIR)/PedCount.pl \
	  $(BINDIR)/kf.pm \
	  $(BINDIR)/Kelvin

.PHONY : kelvin
kelvin : kelvin-$(VERSION)

kelvin-$(VERSION) : libs $(KOBJS) $(OBJS) $(INCS)
	$(CC) -o $@ $(KOBJS) $(OBJS) -lped -lconfig -lklvnutls -lm -lpthread $(LDFLAGS) $(CFLAGS) $(EXTRAFLAG)
	cp $@ $@-$(SVNVERSION)


kelvin.platform : clean check_dist_flag kelvin
	cp kelvin-$(VERSION) bin/kelvin.$(PLATFORM)

check_dist_flag :
ifeq ("$(findstring -DDISTRIBUTION,$(CFLAGS))", "")
	echo You must enable the DISTRIBUTION compilation flag for prebuilt executables!
	exit 1
endif


.PHONY : seq_update/calc_updated_ppl
seq_update/calc_updated_ppl :
	+make -C seq_update -f Makefile calc_updated_ppl
	cp $@ bin/calc_updated_ppl.$(PLATFORM)

%.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) -DVERSION='"V$(VERSION)"' -DSVNVERSION='"$(SVNVERSION)"' $< -o $@

.PHONY : libs
libs :
	+make -C utils -f Makefile all
	+make -C config -f Makefile all
	+make -C pedlib -f Makefile all
ifeq ($(strip $(USE_STUDYDB)), yes)
	+make -C database -f Makefile all
endif

.PHONY : clean
clean :
	make -C pedlib -f Makefile clean
	make -C config -f Makefile clean
	make -C utils -f Makefile clean
ifeq ($(strip $(USE_STUDYDB)), yes)
	make -C database -f Makefile clean
endif
	make -C seq_update -f Makefile clean
	rm -f $(KOBJS) $(OBJS) kelvin-$(VERSION) seq_update/calc_updated_ppl lib/libconfig.a lib/klvnutls.a lib/libped.a
ifeq ($(strip $(USE_STUDYDB)), yes)
	rm -f lib/klvndb.a
endif
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

$(BINDIR)/kelvin.$(PLATFORM) :
ifeq (,$(wildcard bin/kelvin.$(PLATFORM)))
	echo Platform-specific prebuilt executable bin/kelvin.$(PLATFORM) does not exist!
	exit 1
else
	install -o $(OWNER) -g $(GROUP) -m 0755 -p bin/kelvin.$(PLATFORM) $(BINDIR)/kelvin.$(PLATFORM)
	install -o $(OWNER) -g $(GROUP) -m 0755 -p bin/kelvin.$(PLATFORM) $(BINDIR)/kelvin-$(VERSION)
	install -o $(OWNER) -g $(GROUP) -m 0755 -p bin/calc_updated_ppl.$(PLATFORM) $(BINDIR)/calc_updated_ppl-$(VERSION)
endif

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
	perl -pe "s|NO_KELVIN_ROOT|$(ABSBINDIR)|; s|NO_KELVIN_BINARY|$(ABSBINDIR)/kelvin-$(VERSION)|; s|NO_SEQUPDATE_BINARY|$(ABSBINDIR)/calc_updated_ppl|; s|NO_PEDCOUNT_SCRIPT|$(ABSBINDIR)/PedCount.pl|" Kelvin > Kelvin.local
	install -o $(OWNER) -g $(GROUP) -m 0755 -p Kelvin.local $(BINDIR)/Kelvin
