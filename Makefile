# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes
# only.

## Directory into which compiled executables and scripts will be installed.
## For LKS, this should generally NOT be a directory in $PATH.
ifndef BINDIR
  BINDIR=/usr/local/share/kelvin
endif

## Directory on your system's $PATH into which Kelvin-LKS front end scripts 
## will be symlinked
PATHDIR=/usr/local/bin

## User and group IDs by which installed executables and scripts will be owned.
OWNER=root
GROUP=root

## The C compiler to be used to build executables. Pick one.
## GCC (GNU C Compiler)
CC := gcc
## ICC (Intel C Compiler)
# CC := icc

## GCC optimization level, 0=none, 1=default, 2=some (recommended), 3=all
GCCOPT := 2

## LKS Portable Edition site-wide configuration
## These values can be changed after installation by editing
## $(BINDIR)/site_settings.sh
# MySQL distribution location, for setting up transient MySQL servers
LKSPE_MYSQL_BASE="/opt/mysql-5.6"

# Directories that need to also go into $TOOLPATH. This is for if you had more
# than one directory in DEFAULT_TOOLPATH in old versions of Portable Edition;
# typically it's only used at sites that also need to include $MYSQL_BASE/bin
LKSPE_EXTENDED_TOOLPATH=""

# qsub command line arguments used to specify that a job will itself submit
# jobs. If left blank, then it's assumed all jobs can submit subsequent jobs.
LKSPE_JOB_SUBMITS_JOBS=$(shell if [[ "$$HOSTNAME" == "Levi-Montalcini"* ]]; then echo "-q headnode -P headnode"; fi)



##                                                     ##
## Should be no need to make changes beyond this point ##
##                                                     ##



## Command used to download files from the Web. Used to build patched Merlin.
## This is necessary because Merlin's licensing terms do not allow us to
## distribute a patched source tree directly.
#WGET := wget
WGET := curl -L -O
## URL path to the Merlin downloadable tarball
MERLIN_URL := http://csg.sph.umich.edu/abecasis/merlin/download
## Filename of the Merlin downloadable tarball
MERLIN_TARBALL := merlin-1.1.2.tar.gz

## Current Merlin patch file
MERLIN_PATCH := $(shell cd merlin; /bin/ls merlin-*.diff.gz)

## Enable OpenMP support. Requires icc or gcc 4.2+, and GSL. By default, it's
## enabled if it's available.
USE_OPENMP ?= $(shell if [ "$$(cpp --version | head -1 | sed -e 's/[^0-9. ]*//g' -e 's/^ *//' -e 's/ .*//')" \< "4.1" ]; then echo 'no'; else echo 'yes'; fi)

## Enable use of GSL (GNU Scientific Library). Don't forget to set
## INCDIR and LIBDIR (above) accordingly.
USE_GSL := yes

## Enable use of ptmalloc3. Don't forget to set LIBDIR (below) accordingly.
## Not available on OSX. By default, it's enabled if it's available.
USE_PTMALLOC3 ?= $(shell if ldconfig -p 2>/dev/null | grep ptmalloc3 >/dev/null; then echo "yes"; else echo "no"; fi)

## Enable use of Hoard. Don't forget to set LIBDIR (below) accordingly.
# USE_HOARD := yes


VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
LKS_VERSION := $(shell echo `cat .maj`.`cat .min`.`cat .pat`-`cat .lks`)
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

# Repository where we can find the pipeline scripts
# Note that this will probably only work internally for obvious reasons!
PIPELINE_SVNROOT:=https://hodgkin/svn/bcmmtools

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


## Debugging and diagnostic flags. For BCMM use only!

## DISTRIBUTION - Eliminates all diagnostics for distribution. Do not change
## this, as it is a dist search target!
#FILE_CFLAGS += -DDISTRIBUTION

## USESHM and BACKTRACE - Shared memory diagnostics. Backtrace may not always
## be supported and so is conditional
#FILE_CFLAGS += -DUSESHM
#ifneq (,$(wildcard /usr/include/execinfo.h))
#  FILE_CFLAGS += -DBACKTRACE
#endif

## MEMSTATUS and MEMGRAPH - Reports information on elapsed time and memory
## consumption every 30 seconds for diagnostic/reporting purposes. MEMSTATUS
## displays this information directly; MEMGRAPH logs it to a file named
## `kelvin_<PID>_memory.dat` that can be graphed by gnuplot with a simple `load
## "<filename>"` command. 
#FILE_CFLAGS += -DMEMSTATUS
#FILE_CFLAGS += -DMEMGRAPH

## DMUSE - experimental in-house static memory management
#FILE_CFLAGS += -DDMUSE

## DMTRACK - enable exhaustive memory management tracking.
#FILE_CFLAGS += -DDMTRACK

## POLYSTATISTICS - Enables a dump of extensive polynomial build statistics at
## every 8 million raw terms and at build milestones.
#FILE_CFLAGS += -DPOLYSTATISTICS

## POLYCHECK_DL - Do a verification comparison between in-memory built
## polynomials and compiled dynamic-library polynomials (keep them both and
## compare results)
#FILE_CFLAGS += -DPOLYCHECK_DL
#FILE_LDFLAGS += -lsocket -lnsl # needed for this under Solaris

## VERIFY_GSL - Do a verification comparison between our own internal
## statistical routines and those from GSL. Routines will return internal-
## routine result and print a warning if error > 1e-13. Incompatible with
## OpenMP as internal routines are not threadsafe.
#FILE_CFLAGS += -DVERIFY_GSL

## TELLRITA and FULLLOG - Saves all log messages. TELLRITA sends them via UDP
## port 4590 to our head node (hardcoded destination); FULLLOG writes them to
## a file named `kelvin.full_log`.
#FILE_CFLAGS += -DTELLRITA
#FILE_CFLAGS += -DFULLLOG


## Enable use of MySQL database for HLOD storage/retrieval/distributed
## processing (the "likelihood server").
# USE_STUDYDB := yes


## Polynomial compilation flags. Used by special Makefile targets; most of the
## time, you'll just want to run `make specialty` to take advantage of these.

## POLYUSE_DL - Use dynamic libraries built for named polynomials if they're
## present in the current default directory. See olddoc/compiledPolys.html for
## usage details.
#FILE_CFLAGS += -DPOLYUSE_DL
#FILE_LDFLAGS += -ldl

## POLYCODE_DL - Generate C code for polynomials that can be compiled into
## dynamic libraries for later reuse by a Kelvin binary with POLYUSE_DL
## enabled. See olddoc/compiledPolys.html for usage details.
#FILE_CFLAGS += -DPOLYCODE_DL
#FILE_LDFLAGS += -ldl

## POLYCOMP_DL - Compile any generated C code for dynamic polynomial libraries
## on the fly. (They can then be used if POLYUSE_DL is turned on.) See 
## olddoc/compiledPolys.html for usage details.
#FILE_CFLAGS += -DPOLYCOMP_DL
#FILE_LDFLAGS += -ldl

## FAKEEVALUATE - Don't evaluate any generated polynomials. Results will be
## wrong! Used with POLYCODE_DL to do code generation for dynamic polynomial
## libraries.
#FILE_CFLAGS += -DFAKEEVALUATE

## USE_SSD - Highly experimental use of solid state drive when building
## polynomials that are larger than available memory. Recommended only for
## generating code for polynomials. NOT THREAD-SAFE!
#FILE_CFLAGS += -DUSE_SSD



ifeq ($(strip $(USE_STUDYDB)), yes)
  FILE_CFLAGS += -DSTUDYDB $(shell mysql_config --include)
  FILE_LDFLAGS += -lklvndb $(shell mysql_config --libs)
endif

# If GSL support has been enabled, in any fashion
ifeq ($(strip $(USE_GSL)), yes)
  FILE_CFLAGS += -DUSE_GSL
  GSL_FLAGS := yes
endif
ifneq (,$(findstring VERIFY_GSL,$(FILE_CFLAGS)))
  GSL_FLAGS := yes
endif
ifeq ($(strip $(GSL_FLAGS)), yes)
  FILE_CFLAGS += $(shell pkg-config --cflags gsl)
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

.SECONDEXPANSION: 
# this is necessary because otherwise references to $(bindir_pipeline) and $(bindir_lks) won't work as prereqs!

all : kelvin-$(VERSION) seq_update/$(CALC_UPDATED_PPL)

specialty : kelvin-$(VERSION)-no_GSL \
	kelvin-$(VERSION)-POLYUSE_DL \
	kelvin-$(VERSION)-POLYCOMP_DL \
	kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD \
	kelvin-study \
	kelvin-normal

dist : all-pipeline-scripts
	- rm -rf kelvin-$(LKS_VERSION)
	mkdir kelvin-$(LKS_VERSION)
	#mkdir kelvin-$(LKS_VERSION)/bin
	#cp -a bin/{kelvin,$(CALC_UPDATED_PPL)}.* kelvin-$(LKS_VERSION)/bin
	cp -a README .maj .min .pat .lks .svnversion Kelvin Kelvin*.pm CHANGES COPYRIGHT convertconfig.pl *.[ch] compileDL.sh kinfo.pl kelvin-$(LKS_VERSION)
	perl -pe "s|#FILE_CFLAGS \+\= \-DDISTRIBUTION|FILE_CFLAGS \+\= \-DDISTRIBUTION|;" Makefile > kelvin-$(LKS_VERSION)/Makefile
	mkdir kelvin-$(LKS_VERSION)/{lib,utils,pedlib,config,seq_update,database,merlin,LKS}
	cp -a utils/Makefile utils/*.{c,h,pl} kelvin-$(LKS_VERSION)/utils
	cp -a pedlib/Makefile pedlib/*.[ch] kelvin-$(LKS_VERSION)/pedlib
	cp -a config/Makefile config/*.[ch] kelvin-$(LKS_VERSION)/config
	cp -a database/Makefile database/*.[ch] kelvin-$(LKS_VERSION)/database
	cp -a seq_update/Makefile seq_update/*.{c,h,pl} kelvin-$(LKS_VERSION)/seq_update
	cp -a merlin/$(MERLIN_PATCH) kelvin-$(LKS_VERSION)/merlin
	cp -a pipeline-scripts kelvin-$(LKS_VERSION)/
	cp -a auxbin kelvin-$(LKS_VERSION)/
	cd LKS/; cp -a $(lks_target_names) ../kelvin-$(LKS_VERSION)/LKS/; cd ..
	mkdir -p kelvin-$(LKS_VERSION)/test-suite/dynamic-grid/PE/SA_DT
	cp -a test-suite/Makefile kelvin-$(LKS_VERSION)/test-suite
	cp -a test-suite/dynamic-grid/Makefile kelvin-$(LKS_VERSION)/test-suite/dynamic-grid
	cp -a test-suite/dynamic-grid/PE/Makefile kelvin-$(LKS_VERSION)/test-suite/dynamic-grid/PE
	cp -a test-suite/dynamic-grid/PE/SA_DT/* kelvin-$(LKS_VERSION)/test-suite/dynamic-grid/PE/SA_DT
	mkdir -p kelvin-$(LKS_VERSION)/test-suite/seq_update/d-2pt-le
	cp -a test-suite/seq_update/Makefile kelvin-$(LKS_VERSION)/test-suite/seq_update
	cp -a test-suite/seq_update/d-2pt-le/* kelvin-$(LKS_VERSION)/test-suite/seq_update/d-2pt-le
	mkdir -p kelvin-$(LKS_VERSION)/test-suite/frontend/KelvinClasses
	cp -a test-suite/frontend/Makefile kelvin-$(LKS_VERSION)/test-suite/frontend
	cp -a test-suite/frontend/KelvinClasses/* kelvin-$(LKS_VERSION)/test-suite/frontend/KelvinClasses
	mkdir kelvin-$(LKS_VERSION)/doc
	doc/convert_docs.py
	cp -a doc/*.md doc/*.html doc/*.jpg kelvin-$(LKS_VERSION)/doc
	#html2text -rcfile doc/html2text.rc doc/index.html >kelvin-$(LKS_VERSION)/doc/index.txt
	#html2text -rcfile doc/html2text.rc doc/Directives.html >kelvin-$(LKS_VERSION)/doc/Directives.txt
	#html2text -rcfile doc/html2text.rc doc/Allocators.html >kelvin-$(LKS_VERSION)/doc/Allocators.txt
	#html2text -rcfile doc/html2text.rc doc/SSDSupport.html >kelvin-$(LKS_VERSION)/doc/SSDSupport.txt
	#html2text -rcfile doc/html2text.rc doc/compiledPolys.html >kelvin-$(LKS_VERSION)/doc/compiledPolys.txt
	#cp -a doc/*.{html,png,gif} kelvin-$(LKS_VERSION)/doc
	tar -cvzf kelvin-$(LKS_VERSION).tar.gz kelvin-$(LKS_VERSION)/
	rm -rf kelvin-$(LKS_VERSION)

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

install-lks : \
	$(BINDIR)/kelvin-$(VERSION) \
	$(BINDIR)/KelvinConfig.pm \
	$(BINDIR)/KelvinDataset.pm \
	$(BINDIR)/KelvinFamily.pm \
	$(BINDIR)/KelvinIO.pm \
	$(BINDIR)/Kelvin \
	$(BINDIR)/kelvin-study \
	$$(bindir_pipeline) \
	$$(bindir_lks) \
	$(BINDIR)/ready_kelvin_lks_analysis.sh \
	$(BINDIR)/site_settings.sh \
	$(BINDIR)/BCMM \
	$(BINDIR)/kinfo.pl \
	$(BINDIR)/wordDiff.pl \
	$(BINDIR)/merlin \
	$(BINDIR)/minx \
	$(BINDIR)/McSample.run \
	$(BINDIR)/makeped
	ln -s $(BINDIR)/Kelvin $(PATHDIR)/Kelvin
# that rule is because we don't want to insert this bit into the existing
# $(BINDIR)/Kelvin rule, as that would break standalone builds

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
	make -C database -f Makefile clean
ifeq ($(strip $(USE_STUDYDB)), yes)
	make -C database -f Makefile clean
endif
	rm -f $(OBJS) kelvin-$(VERSION) kelvin-$(VERSION)-$(SVNVERSION) seq_update/$(CALC_UPDATED_PPL) lib/libconfig.a lib/klvnutls.a lib/libped.a
ifeq ($(strip $(USE_STUDYDB)), yes)
	rm -f lib/klvndb.a
endif

.PHONY : clean
clean : clean-build clean-merlin clean-pscripts
	make -C seq_update -f Makefile clean
	make -C test-suite -f Makefile clean
	rm -f kelvin-$(VERSION)-no_GSL kelvin-$(VERSION)-POLYUSE_DL kelvin-$(VERSION)-POLYCOMP_DL kelvin-$(VERSION)-POLYCODE_DL_FAKEEVALUATE_SSD kelvin-$(VERSION)-study kelvin-study kelvin-$(VERSION)-normal kelvin-normal
	rm -f LKSPortableEdition-*.tar.gz

.PHONY : clean-pscripts
clean-pscripts : 
ifeq ("$(findstring -DDISTRIBUTION,$(CFLAGS))", "")
	rm -rf pipeline-scripts
endif

.PHONY : clean-merlin
clean-merlin :
	mv merlin merlin-cleaning
	mkdir merlin
	mv merlin-cleaning/$(MERLIN_PATCH) merlin
	rm -rf merlin-cleaning

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


# Kelvin-LKS Portable Edition BCMMTools Pipeline Scripts

# These are all scripts that are currently maintained in a SVN repository other
# than the one the Kelvin tree is in. Therefore, we need to retrieve them.
# 
# The idea here is that we export all our required scripts to a
# pipeline-scripts subdirectory - that's how we "make" each target script.
# Then, when it's time to install, we copy from that directory into $BINDIR.
# Conceptually simple, but surprisingly awkward to implement in make.
# 
# It also doesn't help that there's two directories in the repository that, for
# Portable Edition purposes, we want to convert into one...

# Pipeline scripts included in SVN bcmmtools/
bt_scriptnames := \
    qtdt_munge.pl \
    mysql_run_script.pl \
    kdtk.pl \
    BCMMTools.pm \
    depcheck_funcs.sh \
    depcheck.pl
# Pipeline scripts included in SVN bcmmtools/cleaning/
btc_scriptnames := \
    cleaning_arrayjob.sh \
    cleaning_common.sh \
    error_handling.sh \
    cmdline_parser.sh \
    template_control_analysis_portable.sh \
    template_settings_portable.sh \
    get_ped_complexity.sh \
    prep_for_kelvin.sh \
    setup_for_LKS.sh \
    MCMC_setup.sh \
    MCMC_sizing.sh \
    MCMC_size_pedigree.sh \
    MCMC_sampling.sh \
    MCMC_sample_pedigree.sh \
    MCMC_gathering.sh \
    collect_MCMC_pedigrees.pl \
    source_a_script.sh \
    transient_MySQL_LKS_setup.sh \
    initial_LKS_runs.sh \
    MCMC_kelvin_MP_only_setup.sh \
    LKS_iteration.sh \
    MCMC_kelvin_server_pedigree.sh \
    LKS_pool_ppl.sh \
    kelvin-split-client.pl \
    LKS_transient_servers.sh \
    transient_server_common_custopts.sh \
    mp_marker_strip.pl \
    mysql-transient-ramdisk.cnf \
    MCMC_size_reporting.pl \
    makeped.pl
# Pipeline scripts included in SVN bcmmtools/cleaning/ that require special
# handling at install time
btc_scriptspecials := \
    BCMM \
    ready_kelvin_lks_analysis.sh \
    site_settings.sh
# There's also a pair of tarballs, each with one of the Acid Tests; those are
# handled below.

# The following are a few mutations of the above lists that we'll need later.
# 
# Convert the lists above to "pipeline-scripts/<scriptname>" lists (our
# "intermediate" target names):
bt_targets = $(patsubst %,pipeline-scripts/%,$(bt_scriptnames))
btc_targets = $(patsubst %,pipeline-scripts/%,$(btc_scriptnames) $(btc_scriptspecials))
# Convert the lists above to a "$(BINDIR)/<scriptname>" list (our installation
# target names):
bindir_pipeline = $(patsubst %,$(wildcard $(BINDIR))/%,$(bt_scriptnames) $(btc_scriptnames) sa_dt-acid-test.tgz merlin-only-sadt-acid-test.tgz)

# How our scripts and/or acid test tarballs are brought in
$(bt_targets) :
	@mkdir -p pipeline-scripts
	svn export $(PIPELINE_SVNROOT)/$(notdir $@) $@

$(btc_targets) :
	@mkdir -p pipeline-scripts
	svn export $(PIPELINE_SVNROOT)/cleaning/$(notdir $@) $@

pipeline-scripts/%.tgz :
	@mkdir -p pipeline-scripts
	svn export $(PIPELINE_SVNROOT)/cleaning/test-suite/depcheck/$*
	tar zcf $@ $*
	rm -rf $*

# How they're installed
# (one generic rule, and three special cases)
# FIXME: The below will fail if $(BINDIR) does not exist!
$(bindir_pipeline): $(wildcard $(BINDIR))/%: pipeline-scripts/% $(BINDIR)
	install -o $(OWNER) -g $(GROUP) -m $(shell if [ -x $< ]; then echo "0755"; else echo "0644"; fi) -p $< $@
# special case #1: ready_kelvin_lks_analysis.sh needs to be symlinked to
# somewhere in $PATH
$(BINDIR)/ready_kelvin_lks_analysis.sh : pipeline-scripts/ready_kelvin_lks_analysis.sh
	install -o $(OWNER) -g $(GROUP) -m 0755 -p $< $@
	ln -s $@ $(PATHDIR)
# special case #2: site_settings.sh actually gets modified on the fly with
# variables defined in this Makefile
$(BINDIR)/site_settings.sh : pipeline-scripts/site_settings.sh
	cat $< | sed -e 's|MYSQL_BASE-".*"|MYSQL_BASE-"$(LKSPE_MYSQL_BASE)"|' -e 's|EXTENDED_TOOLPATH-".*"|EXTENDED_TOOLPATH-"$(LKSPE_EXTENDED_TOOLPATH)"|' -e 's|JOB_SUBMITS_JOBS-".*"|JOB_SUBMITS_JOBS-"$(LKSPE_JOB_SUBMITS_JOBS)"|' > sitesettings.local
	install -o $(OWNER) -g $(GROUP) -m 0755 -p sitesettings.local $@
	rm sitesettings.local
# special case #3: BCMM is actually a directory; we need to install that
# directory's contents
$(BINDIR)/BCMM : pipeline-scripts/BCMM
	mkdir -p $@
	install -o $(OWNER) -g $(GROUP) -m 0755 -p $</CLIParser.pm $@/CLIParser.pm

# Generic target for all of the above
all-pipeline-scripts : $(bt_targets) $(btc_targets) pipeline-scripts/sa_dt-acid-test.tgz pipeline-scripts/merlin-only-sadt-acid-test.tgz
	mkdir -p pipeline-scripts
# and a target to create the directory
# Got rid of this because no matter what I tried, make would perpetually redo
# everything in the directory if I tried to make this a dependency. :(
#pipeline-scripts/.dirstamp :
#	mkdir pipeline-scripts && touch $@

# Kelvin-LKS Portable Edition LKS scripts from this repository

# Pipeline scripts included in LKS/ in this repository
# We define a list here because LKS/ has a bunch of other related bits in it
# that aren't needed for Portable Edition but are still necessary for other
# circumstances.
lks_target_names := \
    LKS_transient_server.sql \
    LKS_transient_database.sql \
    LKS_setup_tables.sql \
    LKS_setup_trigger_proc.sql \
    DModelParts.sql \
    QModelParts.sql \
    InitStudy.pl \
    remove_study.sh \
    merlin_lk_prepare.pl \
    merlin_lk_server.pl
bindir_lks = $(patsubst %,$(BINDIR)/%, $(lks_target_names))
$(bindir_lks) : $(wildcard $(BINDIR))/%: LKS/%
	install -o $(OWNER) -g $(GROUP) -m $(shell if [ -x $< ]; then echo "0755"; else echo "0644"; fi) -p $< $@
# Other such scripts located in different spots
$(BINDIR)/wordDiff.pl :
	install -o $(OWNER) -g $(GROUP) -m 0755 -p utils/wordDiff.pl $(BINDIR)/wordDiff.pl
$(BINDIR)/kinfo.pl :
	install -o $(OWNER) -g $(GROUP) -m 0755 -p kinfo.pl $(BINDIR)/kinfo.pl


# Targets for patched Merlin
# This is also needed for Kelvin-LKS Portable Edition.
$(BINDIR)/merlin : merlin/executables/merlin
	install -o $(OWNER) -g $(GROUP) -m 0755 -p merlin/executables/merlin $(BINDIR)/merlin

$(BINDIR)/minx : merlin/executables/minx
	install -o $(OWNER) -g $(GROUP) -m 0755 -p merlin/executables/minx $(BINDIR)/minx

$(BINDIR)/pedstats : merlin/executables/pedstats
	install -o $(OWNER) -g $(GROUP) -m 0755 -p merlin/executables/pedstats $(BINDIR)/pedstats

merlin/executables/% : merlin/Makefile
	+make -C merlin -f Makefile BINDIR=executables executables/$*

# To get the Makefile, we download the tarball from the Merlin distribution
# website, untar it, and patch the result. This is necessary because Merlin's
# licensing terms nominally forbid redistribution of any sort.
merlin/Makefile :
	cd merlin; \
	$(WGET) $(MERLIN_URL)/$(MERLIN_TARBALL) ; \
	tar zxf $(MERLIN_TARBALL) --strip=1 ; \
	tar zxf pedstats*.tar.gz ; \
	zcat $(MERLIN_PATCH) | patch -p1


# Auxillary binaries
# Will work out a better way of distributing these later.
$(BINDIR)/makeped :
	install -o $(OWNER) -g $(GROUP) -m 0755 -p auxbin/makeped $(BINDIR)/makeped

$(BINDIR)/McSample.run :
	install -o $(OWNER) -g $(GROUP) -m 0755 -p auxbin/McSample.run $(BINDIR)/McSample.run

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

uninstall-lks :
	rm -rf $(BINDIR)/kelvin-$(VERSION) \
	$(BINDIR)/KelvinConfig.pm \
	$(BINDIR)/KelvinDataset.pm \
	$(BINDIR)/KelvinFamily.pm \
	$(BINDIR)/KelvinIO.pm \
	$(BINDIR)/kelvin-study \
	$$(bindir_pipeline) \
	$$(bindir_lks) \
	$(BINDIR)/ready_kelvin_lks_analysis.sh \
	$(BINDIR)/site_settings.sh \
	$(BINDIR)/BCMM \
	$(BINDIR)/kinfo.pl \
	$(BINDIR)/wordDiff.pl \
	$(BINDIR)/kinfo.pl \
	$(BINDIR)/merlin \
	$(BINDIR)/minx \
	$(BINDIR)/McSample.run \
	$(BINDIR)/makeped \
	$(BINDIR)/Kelvin
