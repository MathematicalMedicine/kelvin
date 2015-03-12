# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes
# only.

## Variables and options in this Makefile are fully documented in
## doc/compileoptions.html

## Directory into which compiled executables and scripts will be installed.
ifndef BINDIR
  BINDIR=~/mykelvin
endif

## User and group IDs by which installed execuatbles and scripts will be owned.
OWNER=root
GROUP=root

## The C compiler to be used to build executables. Pick one.
## GCC (GNU C Compiler)
CC := gcc
## ICC (Intel C Compiler)
# CC := icc

## GCC optimization level, 0=none, 1=default, 2=some (recommended), 3=all
GCCOPT := 2

##                                                     ##
## Should be no need to make changes beyond this point ##
##                                                     ##

## Enable OpenMP support. Requires icc or gcc 4.2+, and GSL
USE_OPENMP := $(shell ./specialty-scripts/use_openmp.sh)

## Enable use of GSL (GNU Scientific Library). Don't forget to set
## INCDIR and LIBDIR (above) accordingly.
USE_GSL := yes

## Enable use of ptmalloc3. Don't forget to set LIBDIR (below) accordingly.
## Not available on OSX.
USE_PTMALLOC3 := $(shell ./specialty-scripts/use_ptmalloc3.sh)

## Enable use of Hoard. Don't forget to set LIBDIR (below) accordingly.
# USE_HOARD := yes

## Enable use of MySQL database for HLOD storage/retrieval/distributed processing.
# USE_STUDYDB := yes

VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
PLATFORM_NAME := $(shell echo `uname -m`-`uname -s`)
empty:=
space:= $(empty) $(empty)
PLATFORM = $(subst $(space),$(empty),$(PLATFORM_NAME))
KVNLIBDIR := $(shell pwd)/lib
KELVIN_ROOT := $(shell pwd)
TEST_KELVIN := $(KELVIN_ROOT)/kelvin-$(VERSION)
KELVIN_SCRIPT := $(KELVIN_ROOT)/Kelvin
# If the build is occurring under cygwin, then we have to change the name of this binary due to Microsoft IDioT (Installer Detection Technology)
ifeq ("$(findstring CYGWIN,$(PLATFORM))", "")
  CALC_UPDATED_PPL := calc_updated_ppl
else
  CALC_UPDATED_PPL := calc_updtd_ppl
endif
SEQUPDATE_BINARY := $(KELVIN_ROOT)/seq_update/$(CALC_UPDATED_PPL)

## Directories in which optional header files and libraries can be found (GSL,
## etc). Remember you can specify these as command-line macros, e.g. at OSC:
## $ make INCDIR=/home/ccri0005/include LIBDIR=/home/ccri0005/lib
ifndef LIBDIR
  LIBDIR=/usr/local/lib
endif
ifndef INCDIR
  INCDIR=/usr/local/lib/include
endif

ABSBINDIR=$(shell echo $(BINDIR))

# FILE_CFLAGS and FILE_LDFLAGS are built internally based upon the configuration options
FILE_CFLAGS := -O$(GCCOPT) -DGCCOPT=$(GCCOPT)
FILE_LDFLAGS := -dynamic -L$(LIBDIR)

## Compiler warnings. In development, abort on warning
# FILE_CFLAGS += -Wshadow -Werror
FILE_CFLAGS += -Wall

## Enable debugging symbols. Inflicts a small drag (10%) on performance
FILE_CFLAGS += -g

# Required to get re-entrant routines under Solaris, benign on other platforms
FILE_CFLAGS += -D_REENTRANT


# If OpenMP support has been enabled, GSL is required. The GSL-replacement 
# routines are not thread-safe.
ifeq ($(strip $(USE_OPENMP)), yes)
  USE_GSL := yes
  ifeq ($(strip $(CC)), gcc)
    # Compiler flags for GCC
    FILE_CFLAGS += -fopenmp
    FILE_LDFLAGS += -fopenmp
  else
    ifeq ($(strip $(CC)), icc)
      # Compiler flags for ICC
      FILE_CFLAGS += -openmp 
      FILE_LDFLAGS += -openmp
    endif
  endif
endif

# If ptmalloc3 support has been enabled
ifeq ($(strip $(USE_PTMALLOC3)), yes)
  FILE_LDFLAGS += -lptmalloc3
endif

# libpthread is required for the timing thread and various options
FILE_LDFLAGS += -lpthread

# If Hoard support has been enabled
ifeq ($(strip $(USE_HOARD)), yes)
  FILE_LDFLAGS += -lhoard
endif

# If we're building in an svn-managed context, get AND preserve the latest svn version
SVNVERSION := $(subst exported,,$(shell svnversion 2>/dev/null))
ifeq ("$(strip $(SVNVERSION))","")
  SVNVERSION := $(shell cat .svnversion)
else
  UPDATE_SVNVERSION := $(shell echo $(SVNVERSION) >.svnversion)
endif

INCFLAGS := -I$(INCDIR)

# cygwin, testmac? and flair don't recognize the -rdynamic bit...
FILE_LDFLAGS += -rdynamic

# Flags for BCMM use only
# Several of these are used for generating the "specialty" binaries.

#FILE_CFLAGS += -DDISTRIBUTION # Eliminates all diagnostics for distribution, don't change, its a dist search target
#FILE_CFLAGS += -DUSESHM # Enables shared memory diagnostic segment
ifneq (,$(wildcard /usr/include/execinfo.h))
#  FILE_CFLAGS += -DBACKTRACE # Add backtrace where supported
endif
#FILE_CFLAGS += -DMEMSTATUS # Display time and memory consumption every 30 seconds
#FILE_CFLAGS += -DMEMGRAPH # Log terse time and memory consumption info to a data file every 30 seconds for graphing
#FILE_CFLAGS += -DPOLYSTATISTICS # Display extensive polynomial statistics every raw 8Mp and at milestones
#FILE_CFLAGS += -DDMUSE # For our own static memory management, not beneficial as yet.
#FILE_CFLAGS += -DDMTRACK # For our own memory tracking
#FILE_CFLAGS += -DTREEEVALUATE # Use evaluateValue of tree instead of evaluatePoly of list.
#FILE_CFLAGS += -DFAKEEVALUATE # Don't evaluate at all - use only for exercise/compilation. Results will be wrong!
#FILE_CFLAGS += -DPOLYUSE_DL # Dynamically load compiled polynomials for in-process use
#FILE_LDFLAGS += -ldl # ditto
#FILE_CFLAGS += -DPOLYCODE_DL # Enable generation of dynamic library code for selected polynomials
#FILE_CFLAGS += -DPOLYCOMP_DL # Enable compilation of dynamic library code for selected polynomials
#FILE_CFLAGS += -DPOLYCHECK_DL # Keep both built and compiled DL polys and compare results (can be noisy!)
#FILE_CFLAGS += -DTELLRITA # Relay all log messages to rita via UDP
#FILE_CFLAGS += -DFULLLOG # Write all log messages to kelvin.full_log if TELLRITA isn't working
#FILE_LDFLAGS += -lsocket -lnsl # ditto for under Solaris
#FILE_CFLAGS += -DUSE_SSD # Experimental use of solid state drive when building polynomials. NOT THREAD-SAFE!
#FILE_CFLAGS += -DVERIFY_GSL # Use both internal and GSL returning internal and printing if error > 1e-13, no OpenMP
#FILE_LDFLAGS += $(shell pkg-config --libs gsl) # UNCOMMENT THIS if you are doing VERIFY_GSL

ifeq ($(strip $(USE_STUDYDB)), yes)
  FILE_CFLAGS += -DSTUDYDB $(shell mysql_config --include)
  FILE_LDFLAGS += -lklvndb $(shell mysql_config --libs)
endif

# If GSL support has been enabled
ifeq ($(strip $(USE_GSL)), yes)
  FILE_CFLAGS += -DUSE_GSL $(shell pkg-config --cflags gsl)
  FILE_LDFLAGS += $(shell pkg-config --libs gsl)
endif

CFLAGS := $(FILE_CFLAGS) $(ENV_CFLAGS) -DVERSION='"V$(VERSION)"' -DSVNVERSION='"$(SVNVERSION)"'

LDFLAGS := -L$(KVNLIBDIR) $(FILE_LDFLAGS) $(ENV_LDFLAGS)

export KVNLIBDIR VERSION CC CFLAGS LDFLAGS INCFLAGS KELVIN_ROOT TEST_KELVIN KELVIN_SCRIPT SEQUPDATE_BINARY


OBJS = kelvin.o dcuhre.o kelvinInit.o kelvinTerm.o iterationSupport.o integrationSupport.o \
	kelvinHandlers.o kelvinWriteFiles.o dkelvinWriteFiles.o \
	ppl.o saveResults.o trackProgress.o \
	summary_result.o tp_result_hash.o

INCS = kelvin.h kelvinGlobals.h kelvinLocals.h kelvinHandlers.h \
	kelvinInit.h kelvinTerm.h \
	iterationGlobals.h iterationLocals.h iterationSupport.h \
	integrationGlobals.h integrationLocals.h integrationSupport.h \
	kelvinWriteFiles.h dkelvinWriteFiles.h \
	ppl.h dcuhre.h saveResults.h summary_result.h trackProgress.h tp_result_hash.h

all : kelvin-$(VERSION) seq_update/$(CALC_UPDATED_PPL)

specialty : kelvin-$(VERSION)-no_GSL \
	kelvin-$(VERSION)-POLYUSE_DL \
	kelvin-$(VERSION)-POLYCOMP_DL \
	kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD \
	kelvin-study \
	kelvin-normal

dist :
	- rm -rf kelvin-$(VERSION)
	mkdir kelvin-$(VERSION)
	mkdir kelvin-$(VERSION)/bin
	cp -a bin/{kelvin,$(CALC_UPDATED_PPL)}.* kelvin-$(VERSION)/bin
	cp -a README .maj .min .pat .svnversion Kelvin Kelvin*.pm CHANGES COPYRIGHT convertconfig.pl *.[ch] compileDL.sh kelvin-$(VERSION)
	perl -pe "s|#FILE_CFLAGS \+\= \-DDISTRIBUTION|FILE_CFLAGS \+\= \-DDISTRIBUTION|;" Makefile > kelvin-$(VERSION)/Makefile
	mkdir kelvin-$(VERSION)/{lib,utils,pedlib,config,seq_update}
	cp -a utils/Makefile utils/*.{c,h,pl} kelvin-$(VERSION)/utils
	cp -a pedlib/Makefile pedlib/*.[ch] kelvin-$(VERSION)/pedlib
	cp -a config/Makefile config/*.[ch] kelvin-$(VERSION)/config
	cp -a database/Makefile database/*.[ch] kelvin-$(VERSION)/database
	cp -a seq_update/Makefile seq_update/*.{c,h,pl} kelvin-$(VERSION)/seq_update
	mkdir -p kelvin-$(VERSION)/test-suite/dynamic-grid/PE/SA_DT
	cp -a test-suite/Makefile kelvin-$(VERSION)/test-suite
	cp -a test-suite/dynamic-grid/Makefile kelvin-$(VERSION)/test-suite/dynamic-grid
	cp -a test-suite/dynamic-grid/PE/Makefile kelvin-$(VERSION)/test-suite/dynamic-grid/PE
	cp -a test-suite/dynamic-grid/PE/SA_DT/* kelvin-$(VERSION)/test-suite/dynamic-grid/PE/SA_DT
	mkdir -p kelvin-$(VERSION)/test-suite/seq_update/d-2pt-le
	cp -a test-suite/seq_update/Makefile kelvin-$(VERSION)/test-suite/seq_update
	cp -a test-suite/seq_update/d-2pt-le/* kelvin-$(VERSION)/test-suite/seq_update/d-2pt-le
	mkdir -p kelvin-$(VERSION)/test-suite/frontend/KelvinClasses
	cp -a test-suite/frontend/Makefile kelvin-$(VERSION)/test-suite/frontend
	cp -a test-suite/frontend/KelvinClasses/* kelvin-$(VERSION)/test-suite/frontend/KelvinClasses
	mkdir kelvin-$(VERSION)/doc
	html2text -rcfile doc/html2text.rc doc/index.html >kelvin-$(VERSION)/doc/index.txt
	html2text -rcfile doc/html2text.rc doc/Directives.html >kelvin-$(VERSION)/doc/Directives.txt
	html2text -rcfile doc/html2text.rc doc/Allocators.html >kelvin-$(VERSION)/doc/Allocators.txt
	html2text -rcfile doc/html2text.rc doc/SSDSupport.html >kelvin-$(VERSION)/doc/SSDSupport.txt
	html2text -rcfile doc/html2text.rc doc/compiledPolys.html >kelvin-$(VERSION)/doc/compiledPolys.txt
	cp -a doc/*.{html,png,gif} kelvin-$(VERSION)/doc
	tar -cvzf kelvin-$(VERSION).tar.gz kelvin-$(VERSION)/
	rm -rf kelvin-$(VERSION)

install : $(BINDIR) \
	$(BINDIR)/kelvin-$(VERSION) \
          $(BINDIR)/$(CALC_UPDATED_PPL) \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh \
	  $(BINDIR)/KelvinConfig.pm \
	  $(BINDIR)/KelvinDataset.pm \
	  $(BINDIR)/KelvinFamily.pm \
	  $(BINDIR)/KelvinIO.pm \
	  $(BINDIR)/Kelvin

install-prebuilt : $(BINDIR)/kelvin.$(PLATFORM) \
          $(BINDIR)/convert_br.pl \
	  $(BINDIR)/compileDL.sh \
	  $(BINDIR)/$(MODS) \
	  $(BINDIR)/KelvinConfig.pm \
	  $(BINDIR)/KelvinDataset.pm \
	  $(BINDIR)/KelvinFamily.pm \
	  $(BINDIR)/KelvinIO.pm \
	  $(BINDIR)/Kelvin

install-specialty : install \
	$(BINDIR)/kelvin-$(VERSION)-no_GSL \
	$(BINDIR)/kelvin-$(VERSION)-POLYUSE_DL \
	$(BINDIR)/kelvin-$(VERSION)-POLYCOMP_DL \
	$(BINDIR)/kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD \
	$(BINDIR)/kelvin-study \
	$(BINDIR)/kelvin-normal

.PHONY : kelvin
kelvin : kelvin-$(VERSION)

$(BINDIR) : 
	mkdir -p $(BINDIR)

kelvin-$(VERSION) : libs $(OBJS) $(INCS)
	$(CC) -o $@ $(OBJS) $(CFLAGS) -lped -lconfig -lklvnutls $(LDFLAGS) -lm
	cp $@ $@-$(SVNVERSION)


kelvin.platform : clean check_dist_flag kelvin
	cp kelvin-$(VERSION) bin/kelvin.$(PLATFORM)

check_dist_flag :
ifeq ("$(findstring -DDISTRIBUTION,$(CFLAGS))", "")
	echo You must enable the DISTRIBUTION compilation flag for prebuilt executables!
	exit 1
endif


.PHONY : seq_update/$(CALC_UPDATED_PPL)
seq_update/$(CALC_UPDATED_PPL) :
	+make -C seq_update -f Makefile calc_updated_ppl
	# The following dance is the fast way of being able to build both the same and differently
	# named final target without modifying a subordinate makefile. Its harmless, ignore it.
	mv seq_update/calc_updated_ppl seq_update/intermediate-name
	mv seq_update/intermediate-name $@
	cp $@ bin/$(CALC_UPDATED_PPL).$(PLATFORM)

.PHONY : libs
libs :
	+make -C utils -f Makefile all
	+make -C config -f Makefile all
	+make -C pedlib -f Makefile all
ifeq ($(strip $(USE_STUDYDB)), yes)
	+make -C database -f Makefile all
endif

.PHONY : clean-build
clean-build : 
	make -C pedlib -f Makefile clean
	make -C config -f Makefile clean
	make -C utils -f Makefile clean
ifeq ($(strip $(USE_STUDYDB)), yes)
	make -C database -f Makefile clean
endif
	rm -f $(OBJS) kelvin-$(VERSION) kelvin-$(VERSION)-$(SVNVERSION) seq_update/$(CALC_UPDATED_PPL) lib/libconfig.a lib/klvnutls.a lib/libped.a
ifeq ($(strip $(USE_STUDYDB)), yes)
	rm -f lib/klvndb.a
endif

.PHONY : clean
clean : clean-build
	make -C seq_update -f Makefile clean
	make -C test-suite -f Makefile clean
	rm -f kelvin-$(VERSION)-no_GSL kelvin-$(VERSION)-POLYUSE_DL kelvin-$(VERSION)-POLYCOMP_DL kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD kelvin-$(VERSION)-study kelvin-study kelvin-$(VERSION)-normal kelvin-normal

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
	install -o $(OWNER) -g $(GROUP) -m 0755 -p bin/$(CALC_UPDATED_PPL).$(PLATFORM) $(BINDIR)/$(CALC_UPDATED_PPL)-$(VERSION)
endif

$(BINDIR)/$(CALC_UPDATED_PPL) : seq_update/$(CALC_UPDATED_PPL)
	install -o $(OWNER) -g $(GROUP) -m 0755 -p seq_update/$(CALC_UPDATED_PPL) $(BINDIR)/$(CALC_UPDATED_PPL)

$(BINDIR)/convert_br.pl : seq_update/convert_br.pl
	install -o $(OWNER) -g $(GROUP) -m 0755 -p seq_update/convert_br.pl $(BINDIR)/convert_br.pl

$(BINDIR)/compileDL.sh : compileDL.sh
	install -o $(OWNER) -g $(GROUP) -m 0755 -p compileDL.sh $(BINDIR)/compileDL.sh

$(BINDIR)/%.pm : %.pm
	install -o $(OWNER) -g $(GROUP) -m 0644 -p $< $@

$(BINDIR)/Kelvin : Kelvin
	perl -pe "s|NO_KELVIN_ROOT|$(ABSBINDIR)|; s|NO_KELVIN_BINARY|$(ABSBINDIR)/kelvin-$(VERSION)|; s|NO_SEQUPDATE_BINARY|$(ABSBINDIR)/$(CALC_UPDATED_PPL)|;" Kelvin > Kelvin.local
	install -o $(OWNER) -g $(GROUP) -m 0755 -p Kelvin.local $(BINDIR)/Kelvin
	rm Kelvin.local


# Specialty binaries follow.
# These represent four different variants of Kelvin usage:
#  * standalone usage (kelvin-normal)
#  * standalone usage without GSL (kelvin-no_GSL)
#  * compiled polynomial usage (kelvin-POLY*)
#  * Likelihood Server usage (kelvin-study)

# GSL-free version. OpenMP also disabled as it requires GSL.
.PHONY : kelvin-no_GSL
kelvin-no_GSL : kelvin-$(VERSION)-no_GSL

kelvin-$(VERSION)-no_GSL :
	make clean-build
	make USE_GSL=no USE_OPENMP=no kelvin
	mv kelvin-$(VERSION) kelvin-$(VERSION)-no_GSL

$(BINDIR)/kelvin-$(VERSION)-no_GSL : kelvin-$(VERSION)-no_GSL
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION)-no_GSL $(BINDIR)/kelvin-$(VERSION)-no_GSL

# Dynamic polynomial version with OpenMP. Use with existing compiled
# polynomials.
# Set OMP_NUM_THREADS=<something big> for best performance after compilation of
# DLs. Doesn't do compilation if DL not found.
.PHONY: kelvin-POLYUSE_DL
kelvin-POLYUSE_DL : kelvin-$(VERSION)-POLYUSE_DL

kelvin-$(VERSION)-POLYUSE_DL :
	make clean-build
	make ENV_CFLAGS=" -DMEMGRAPH -DPOLYSTATISTICS -DPOLYUSE_DL" ENV_LDFLAGS=" -ldl" kelvin
	mv kelvin-$(VERSION) kelvin-$(VERSION)-POLYUSE_DL

$(BINDIR)/kelvin-$(VERSION)-POLYUSE_DL : kelvin-$(VERSION)-POLYUSE_DL
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION)-POLYUSE_DL $(BINDIR)/kelvin-$(VERSION)-POLYUSE_DL

# Polynomial compiling version - build, code and then compile and evaluate DLs.
# Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
.PHONY: kelvin-POLYCOMP_DL
kelvin-POLYCOMP_DL : kelvin-$(VERSION)-POLYCOMP_DL

kelvin-$(VERSION)-POLYCOMP_DL : 
	make clean-build
	make ENV_CFLAGS=" -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DPOLYCOMP_DL" ENV_LDFLAGS=" -ldl" kelvin
	mv kelvin-$(VERSION) kelvin-$(VERSION)-POLYCOMP_DL

$(BINDIR)/kelvin-$(VERSION)-POLYCOMP_DL : kelvin-$(VERSION)-POLYCOMP_DL
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION)-POLYCOMP_DL $(BINDIR)/kelvin-$(VERSION)-POLYCOMP_DL

# Polynomial generating version - build and code DLs, pretend to evaluate. Best
# for making hay out of limited memory.
# Again, notice no OpenMP, as it gains no speed and loses us lots of memory.
.PHONY: kelvin-POLYCODE_DL_FAKEEVALUATE_SSD
kelvin-POLYCODE_DL_FAKEEVALUATE_SSD : kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD

kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD :
	make clean-build
	make ENV_CFLAGS=" -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DUSE_SSD -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE" ENV_LDFLAGS=" -ldl" kelvin
	mv kelvin-$(VERSION) kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD

$(BINDIR)/kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD : kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD $(BINDIR)/kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD

# Likelihood server version. Requires mySQL version 5 or better!
kelvin-study: kelvin-$(VERSION)-study
	cp kelvin-$(VERSION)-study kelvin-study

kelvin-$(VERSION)-study :
	make clean-build USE_STUDYDB=yes
	make USE_STUDYDB=yes kelvin
	mv kelvin-$(VERSION) kelvin-$(VERSION)-study

$(BINDIR)/kelvin-study : kelvin-study
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION)-study $(BINDIR)/kelvin-$(VERSION)-study
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-study $(BINDIR)/kelvin-study

# Standard standalone version
kelvin-normal : kelvin-$(VERSION)-normal
	cp kelvin-$(VERSION)-normal kelvin-normal

kelvin-$(VERSION)-normal :
	make clean-build
	make kelvin
	mv kelvin-$(VERSION) kelvin-$(VERSION)-normal

$(BINDIR)/kelvin-normal : kelvin-normal
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-normal $(BINDIR)/kelvin-normal
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kelvin-$(VERSION)-normal $(BINDIR)/kelvin-$(VERSION)-normal


# Uninstallation routines follow.

uninstall :
	rm -f $(BINDIR)/kelvin-$(VERSION) \ 
	$(BINDIR)/$(CALC_UPDATED_PPL) \ 
	$(BINDIR)/convert_br.pl \ 
	$(BINDIR)/compileDL.sh \ 
	$(BINDIR)/KelvinConfig.pm \ 
	$(BINDIR)/KelvinDataset.pm \ 
	$(BINDIR)/KelvinFamily.pm \ 
	$(BINDIR)/KelvinIO.pm \ 
	$(BINDIR)/Kelvin

uninstall-prebuilt : 
	rm -f $(BINDIR)/kelvin.$(PLATFORM) \ 
	$(BINDIR)/kelvin-$(VERSION) \ 
	$(BINDIR)/$(CALC_UPDATED_PPL)-$(VERSION) \ 
	$(BINDIR)/convert_br.pl \ 
	$(BINDIR)/compileDL.sh \ 
	$(BINDIR)/$(MODS) \ 
	$(BINDIR)/KelvinConfig.pm \ 
	$(BINDIR)/KelvinDataset.pm \ 
	$(BINDIR)/KelvinFamily.pm \ 
	$(BINDIR)/KelvinIO.pm \ 
	$(BINDIR)/Kelvin

uninstall-specialty :
	rm -f $(BINDIR)/kelvin-$(VERSION)-no_GSL \ 
	$(BINDIR)/kelvin-$(VERSION)-POLYUSE_DL \ 
	$(BINDIR)/kelvin-$(VERSION)-POLYCOMP_DL \ 
	$(BINDIR)/kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD \ 
	$(BINDIR)/kelvin-$(VERSION)-study \ 
	$(BINDIR)/kelvin-study \ 
	$(BINDIR)/kelvin-$(VERSION)-normal \ 
	$(BINDIR)/kelvin-normal
