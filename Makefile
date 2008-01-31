# Multiprocessor Linkage Analysis
# Alberto Maria Segre, Yungui Huang, Nathan Burnette, others?
#
# Copyright (c) 2006, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

CC = gcc

# Base location for NICE installation.
#NICEDIR = /usr/local/nice
NICEDIR = /usr/local

# Compile options:
# -DDEBUG
#     Turns on debugging information.
#
#CFLGS	+= -DNOFLOAT 
CFLGS	+= -Wall			# Always leave this one on
CFLGS	+= -O3	 			# Production level
#CFLGS  += -g  -DDEBUG		# Debug level
#CFLGS   += -DNO_POLYNOMIAL
CFLGS += -pg                   # profile - debugging 
#CFLGS += -pedantic   #memory leak tracing

.EXPORT_ALL_VARIABLES:
######################################################################
# No user changes below this point!
ARCH 	:= $(shell uname -m)
OS  	:= $(shell uname -s)	# Cygwin needs -o?
VERSION	:= $(shell echo `cat .maj`.`cat .min`.`cat .pat`)

# Nice library and header files
NLIBDIR		:= $(NICEDIR)/lib
NINCDIR		:= $(NICEDIR)/include
NBINDIR		:= $(NICEDIR)/bin

# Destination for utilities and libraries
INCDIR := $(PWD)/include
LIBDIR := $(PWD)/lib

CFLGS  += -Wall -I$(INCDIR) -L$(LIBDIR) -I$(NINCDIR) -L$(NLIBDIR)
######################################################################
# File sets. 
NBIN	= kelvin
LIB	= -lped -lutils -lgsl -lgslcblas -lsw -lm 
NLIB	= # -lniceapi -lnicecom -lniceaux
SRC	= kelvin.c config.c ppl.c
INC	= kelvin.h
TAR	= Makefile $(SRC) $(INC) config.c utils pedlib .maj .min .pat .dat doc kelvin.conf

######################################################################
# Determine application name, which depends on version and
# architecture information.
NICEAPP = $(NBIN)-$(VERSION)-$(ARCH)
ifneq ($(NICEOS),Cygwin)
  NICEEXE = $(NBINDIR)/$(NICEAPP)
else
  NICEEXE = $(NBINDIR)/$(NICEAPP).exe
endif

######################################################################
# Make an object file from source
%.o      : %.c $(INCS)
	$(CC) -c $(CFLGS) -o $@ $<

######################################################################
# Manual section.
MSEC	= 1

######################################################################
# Kelvin.
kelvin: Makefile libutils.a libped.a # man
	$(CC) $(CFLGS) $(SRC) -o $@ $(NLIB) $(LIB)

######################################################################
# Default target.  Installs application in $NBINDIR.
.PHONY: install
install: $(NBIN)

######################################################################
# Kelvin utils library.
libutils.a: Makefile
	-@if [ ! -x $(INCDIR) ] ; then mkdir $(INCDIR); fi
	-@if [ ! -x $(LIBDIR) ] ; then mkdir $(LIBDIR); fi
	-@cd utils && $(MAKE) $(MAKECMDGOALS)

######################################################################
# Kelvin pedigree library.
libped.a: Makefile
	-@if [ ! -x $(INCDIR) ] ; then mkdir $(INCDIR); fi
	-@if [ ! -x $(LIBDIR) ] ; then mkdir $(LIBDIR); fi
	-@cd pedlib && $(MAKE) $(MAKECMDGOALS)

######################################################################
.PHONY: man
man: kelvin.man
	@sed -e 's/%%DATE%%/$(DATE)/' < kelvin.man > kelvin.$(MSEC)
#	@gtbl kelvin.$(MSEC) | groff -man -Tps > Manual.ps
#	@ps2pdf Manual.ps
#	@rm Manual.ps

######################################################################
# Installs application in $NBINDIR.
.PHONY: install
install: kelvin #man
	@chmod 755 $(NBIN)
	/bin/cp kelvin $(NBINDIR)/kelvin-$(VERSION)
#	@cp $(NBIN) $(NICEEXE)
#	@-cp kelvin.$(MSEC) /usr/local/man/man$(MSEC)

######################################################################
# Remove installed version; leaves local binary.
.PHONY: remove
remove:
	-@rm $(NICEEXE)

######################################################################
# Utilities.

# To ensure a clean build, use "make fresh."
.PHONY: fresh
fresh: clean $(NBIN)

# Deletes all compiled files.
.PHONY: clean
clean:
	-@rm -f $(NBIN) TAGS *~ $(NBIN)-$(VERSION).{ps,tgz}
	-@cd utils && $(MAKE) $(MAKECMDGOALS)
	-@cd pedlib && $(MAKE) $(MAKECMDGOALS)

.PHONY: ps
ps:
	@enscript -Ec -2r -o $(NBIN)-$(VERSION).ps Makefile kelvin.h kelvin.c config.c
#	@enscript -Ec -2r -o $(NBIN)-$(VERSION).ps Makefile $(INC) $(SRC) 

.PHONY: tags
tags:
	@etags *.[ch]

.PHONY: tgz
tgz: clean
	@tar -czf $(NBIN)-$(VERSION).tgz -C .. $(foreach ENTRY,$(TAR),$(NBIN)-$(VERSION)/$(ENTRY))
	-@ls -lsag $(NBIN)-$(VERSION).tgz

######################################################################
# Creates a new release.
.PHONY: major
major:
	@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NBIN).1
	@cat .maj > /tmp/$(NBIN).2
	@echo "1 + p" >> /tmp/$(NBIN).2
	@dc < /tmp/$(NBIN).2 > .maj
	@echo "0" > .min
	@echo "0" > .pat
	@date > .dat
	@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NBIN).2
	@echo "Version" `cat /tmp/$(NBIN).2`";" `cat .dat`
	@mv ../$(NBIN)-`cat /tmp/$(NBIN).1` ../$(NBIN)-`cat /tmp/$(NBIN).2`
	@cd ../$(NBIN)-`cat /tmp/$(NBIN).2`
	@rm /tmp/$(NBIN).1 /tmp/$(NBIN).2

# Creates a new minor version.
.PHONY: minor
minor:
	@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NBIN).1
	@cat .min > /tmp/$(NBIN).2
	@echo "1 + p" >> /tmp/$(NBIN).2
	@dc < /tmp/$(NBIN).2 > .min
	@echo "0" > .pat
	@date > .dat
	@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NBIN).2
	@echo "Version" `cat /tmp/$(NBIN).2`";" `cat .dat`
	@mv ../$(NBIN)-`cat /tmp/$(NBIN).1` ../$(NBIN)-`cat /tmp/$(NBIN).2`
	@cd ../$(NBIN)-`cat /tmp/$(NBIN).2`
	@rm /tmp/$(NBIN).1 /tmp/$(NBIN).2

# Creates a new patch version.
.PHONY: patch
patch:
	@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NBIN).1
	@cat .pat > /tmp/$(NBIN).2
	@echo "1 + p" >> /tmp/$(NBIN).2
	@dc < /tmp/$(NBIN).2 > .pat
	@date > .dat
	@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NBIN).2
	@echo "Version" `cat /tmp/$(NBIN).2`";" `cat .dat`
	@mv ../$(NBIN)-`cat /tmp/$(NBIN).1` ../$(NBIN)-`cat /tmp/$(NBIN).2`
	@cd ../$(NBIN)-`cat /tmp/$(NBIN).2`
	@rm /tmp/$(NBIN).1 /tmp/$(NBIN).2

# Check variable values.
.PHONY: check
check:
	@echo "$(NBIN) version $(VERSION) for $(OS)"
	@echo NICEDIR=$(NICEDIR)
	@echo NBINDIR=$(NBINDIR)
	@echo INCDIR=$(INCDIR)
	@echo LIBDIR=$(LIBDIR)

