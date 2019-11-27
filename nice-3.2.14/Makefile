# NICE System
# Alberto Maria Segre
#
# Copyright (c) 1997-2005, The University of Iowa.  All rights
# reserved.  Permission is hereby given to use and reproduce this
# software for non-profit educational purposes only.

.EXPORT_ALL_VARIABLES:
######################################################################
# Install NICE system. 
#
# Must be root in order to install "debian" "redhat" or "slackware",
# which make appropriate use of SysV init. When installing "debug" the
# installation will reside directly in the distribution directory;
# you'll have to start your own daemons manually. 
#SYSTYPE	:= debug
#SYSTYPE	:= debian
SYSTYPE	:= redhat
#SYSTYPE	:= slackware
#SYSTYPE	:= cygwin

# Name of distribution directory. Some versions of make (e.g., old
# RedHat distributions and Slackware's for starters) don't expand the
# shell variable properly, so it may be necessary to give the name of
# the directory explicitly.
DISTDIR := $(PWD)

# Name of NICE user/group account, which will be created if not
# already available. NICE daemon and applications run with this user's
# permissions. The user account's password will be locked and no shell
# will be provided.
NICEUSR := niced
NICEGRP := niced
NICEID  := 88

######################################################################
# No user changes below this point!
SHELL		= /bin/sh

NICE		:= nice
NICEVERSION	:= $(shell echo `cat .maj`.`cat .min`.`cat .pat`)
NICEDATE 	:= $(shell echo `cat .dat`)

# Determine the architecture: should typically be i686 or x86_64 for
# Athlon 64 machines. We'd like to use HOSTTYPE but that's a Bash
# shell builtin, which may not be that portable: and while uname -p
# seems a better bet, it sometimes is not set properly, returning
# "unknown."
#
# Note that we want to ensure there are not / or - characters in the 
# architecture name so that, e.g., it isn't misinterpreted as a directory.
NICEARCH 	:= $(shell uname -m | sed -e 's/[\/,-]/+/g')
NICEOS  	:= $(shell uname -s)	# Cygwin needs -o?

######################################################################
# Directories. 
SRCDIR	:= $(DISTDIR)

# Installation prefix. The default is to install in /usr/local, unless
# you are debugging or installing on a Windows machine. Note order of
# the clauses here is important, in case you are in fact debugging on
# a Windows machine!
PREFIX := /usr/local
ifeq ($(SYSTYPE),cygwin)
  PREFIX = "C:/Program Files/Nice $(NICEVERSION)"
endif
ifeq ($(SYSTYPE),debug)
  PREFIX = $(DISTDIR)/local
endif

# Default values for NICEDIR, SBINDIR, BINDIR and CFGDIR
NICEDIR := $(PREFIX)/nice
SBINDIR	:= $(PREFIX)/sbin
BINDIR	:= $(PREFIX)/bin
CFGDIR	:= $(PREFIX)/etc
# Compute values for INCDIR, LIBDIR and NBINDIR
INCDIR	:= $(NICEDIR)/include
ifeq ($(NICEARCH),x86_64)
  LIBDIR0	:= $(NICEDIR)/lib64
else
  LIBDIR0	:= $(NICEDIR)/lib32
endif
LIBDIR	:= $(NICEDIR)/lib
NCFGDIR	:= $(NICEDIR)/etc
NBINDIR	:= $(NICEDIR)/bin

######################################################################
# Make options. Note that the -C must always be last!
# Deprecated: use explicit cd's rather than -C for HP's Make sake.
# The -n flag simply causes make to echo what it would do, without doing it.
#MFLGS	:= -C
#MFLGS	:= -n -C 

# Compiler options. Careful: we're fetching the libraries from the source
# directory within NICE. Applications will presumably fetch them from the 
# NICE installation directory instead.
CC	:= gcc
C++	:= g++
G77	:= g77
CFLGS	:= -Wall -I$(INCDIR) -L$(LIBDIR)

# NODYNAMIC disables dynamic daemon reorganization in niced.
#CFLGS	+= -DNODYNAMIC

# NOALARM disables the software timers, which may interfere with
# debugging. TMRDBG turns on alarm debugging.
#CFLGS	+= -DNOALARM
#CFLGS	+= -DTMRDBG

# NODEFER disables the graceful deferral of incoming messages. 
# Requires ALARMs.
#CFLGS	+= -DNODEFER

# BOUNDED enables a search application's node expansion limits.
#CFLGS	+= -DBOUNDED

# NOSTREAM disables the use of sockread/sockwrite streaming primitives
# in NICE applications.
#CFLGS	+= -DNOSTREAM

# DEBUG turns on application-level debugging, while LIBDBG turns on
# library-level debugging. 
ifeq ($(SYSTYPE),debug)
  CFLGS	+= -g -DDEBUG #-DLIBDBG
else
  CFLGS += -O4
endif

# Ugh! HPUX doesn't handle accept() correctly using the standard CC
# without some additional defines. These would not be necessary if our
# GCC installation on HPUX weren't so irretrievably broken. A similar
# situation holds with the rmdir flags, which would normally just be
# placed elsewhere.
ifeq ($(NICEARCH), 9000+785)
  CC	:= cc
  CFLGS += -D_XOPEN_SOURCE=1 -D_XOPEN_SOURCE_EXTENDED=1 -Ae
  RDFLGS =
else
  RDFLGS = --ignore-fail-on-non-empty
endif

# NOTE: I (Jo) have disabled the below targets, as they all involve
# installation and that absolutely should not be the default.
#######################################################################
## Default target. Install in appropriate locations.
#.PHONY: install
#install: niced ntest account
#	-@if [ ! -f $(CFGDIR)/niced.conf ] ; then cp $(SRCDIR)/etc/niced.conf $(CFGDIR)/niced.conf; echo -e "\nConfigure NICE in $(CFGDIR)/niced.conf\n"; fi
#ifeq ($(SYSTYPE),debian)
#	-@sed -e's%SBINDIR%$(SBINDIR)%g' < $(SRCDIR)/etc/init.debian | sed -e's%NICEUSR%$(NICEUSR)%g' | sed -e's%NICEGRP%$(NICEGRP)%g' > /etc/init.d/niced
#	-@chmod o+x /etc/init.d/niced
#	-@/usr/sbin/update-rc.d niced defaults 90 02
#endif
#ifeq ($(SYSTYPE),redhat)
#	-@sed -e's%SBINDIR%$(SBINDIR)%g' < $(SRCDIR)/etc/init.redhat | sed -e's%NICEUSR%$(NICEUSR)%g' > /etc/init.d/niced
#	-@chmod o+x /etc/init.d/niced
#	-@ln -s /etc/init.d/niced /etc/rc0.d/K02niced
#	-@ln -s /etc/init.d/niced /etc/rc1.d/K02niced
#	-@ln -s /etc/init.d/niced /etc/rc3.d/S90niced
#	-@ln -s /etc/init.d/niced /etc/rc5.d/S90niced
#	-@ln -s /etc/init.d/niced /etc/rc6.d/K02niced
#endif
#ifeq ($(SYSTYPE),slackware)
#	-@sed -e's%SBINDIR%$(SBINDIR)%g' < $(SRCDIR)/etc/init.slack | sed -e's%NICEUSR%$(NICEUSR)%g' > /etc/rc.d/rc.niced
#	-@chmod o+x /etc/rc.d/rc.niced
#	-@if ! grep "niced" /etc/rc.d/rc.local ; then \
#	echo "# Start NICE daemon." >> /etc/rc.d/rc.local ; \
#	echo "if [ -x /etc/rc.d/rc.niced ] ; then" >> /etc/rc.d/rc.local ; \
#	echo "  . /etc/rc.d/rc.niced" >> /etc/rc.d/rc.local ; \
#	echo "fi" >> /etc/rc.d/rc.local ; fi
#endif
#
#######################################################################
## User account. Not needed for debug mode installations since they run
## with the debugging user's privileges.
#.PHONY: account
#ifeq ($(SYSTYPE),debug)
#account: dirs
#	echo "No account created."
#else
#account: dirs
#	-@if ! grep $(NICEGRP) /etc/group > /dev/null ; then /usr/sbin/groupadd -g $(NICEID) $(NICEGRP) ; fi
#	-@if ! grep $(NICEUSR) /etc/passwd > /dev/null ; then /usr/sbin/useradd -u $(NICEID) -s /bin/false -c "NICE daemon" -g $(NICEGRP) -d $(NBINDIR) $(NICEUSR) ; fi
#	-@/usr/sbin/usermod -L -s /bin/false -c "NICE Daemon" -g $(NICEGRP) -d $(NBINDIR) $(NICEUSR)
#	-@/bin/chown -Rf $(NICEUSR).$(NICEGRP) $(NBINDIR)
#endif
#
#######################################################################
## Directories.
#.PHONY: dirs
#dirs: Makefile
#	-@if [ ! -x $(PREFIX) ] ; then mkdir $(PREFIX); fi
#	-@if [ ! -x $(NICEDIR) ] ; then mkdir $(NICEDIR); fi
#	-@if [ ! -x $(INCDIR) ] ; then mkdir $(INCDIR); fi
#	-@if [ ! -x $(LIBDIR0) ] ; then mkdir $(LIBDIR0); fi
#	-@rm $(LIBDIR)
#	-@ln -s $(LIBDIR0) $(LIBDIR)
#	-@if [ ! -x $(BINDIR) ] ; then mkdir $(BINDIR); fi
#	-@if [ ! -x $(SBINDIR) ] ; then mkdir $(SBINDIR); fi
#	-@if [ ! -x $(CFGDIR) ] ; then mkdir $(CFGDIR); fi
#	-@if [ ! -x $(NBINDIR) ] ; then mkdir $(NBINDIR); fi
#
#######################################################################
## Include files.
#.PHONY: inc
#inc: Makefile dirs
#	-@cp $(SRCDIR)/nicelib/niceaux.h $(INCDIR)
#	-@cp $(SRCDIR)/nicelib/nicecom.h $(INCDIR)
#	-@cp $(SRCDIR)/nicelib/niceapi.h $(INCDIR)
#
######################################################################
# Libraries. Depends implicitly on default NICEPORT from config.h.
.PHONY: nicelib config.h
#nicelib: inc $(SRCDIR)/nicelib/Makefile
nicelib: $(SRCDIR)/nicelib/Makefile
	-@cd $(SRCDIR)/nicelib && $(MAKE) $(MAKECMDGOALS)
#	-@make $(MAKECMDGOALS) $(MFLGS) $(SRCDIR)/nicelib

######################################################################
# NICE daemon and control programs. 
.PHONY: niced
niced: nicelib
	-@cd $(SRCDIR)/niced && $(MAKE) $(MAKECMDGOALS)
#	-@make $(MAKECMDGOALS) $(MFLGS) $(SRCDIR)/niced

######################################################################
# NICE-aware test application. 
.PHONY: ntest
ntest: nicelib
	-@cd $(SRCDIR)/ntest && $(MAKE) $(MAKECMDGOALS)
#	-@make $(MAKECMDGOALS) $(MFLGS) $(SRCDIR)/ntest

######################################################################
# Utilities.

# To ensure a clean build, use "make fresh."
.PHONY: fresh
fresh: clean all

# Deletes all compiled files.
.PHONY: clean
clean:
	-@rm -f *~ nice-*.ps nice-*.html *-$(NICEVERSION).tgz *-$(NICEVERSION).tar *-$(NICEVERSION).zip $(SRCDIR)/include/*~
	-@cd $(SRCDIR)/nicelib && $(MAKE) clean
	-@cd $(SRCDIR)/niced && $(MAKE) clean
	-@cd $(SRCDIR)/ntest && $(MAKE) clean
#	-@make $(MFLGS) $(SRCDIR)/nicelib clean
#	-@make $(MFLGS) $(SRCDIR)/niced clean
#	-@make $(MFLGS) $(SRCDIR)/ntest clean

.PHONY: tags
tags:
	-@cd $(SRCDIR)/nicelib && $(MAKE) tags
	-@cd $(SRCDIR)/niced && $(MAKE) tags
	-@cd $(SRCDIR)/ntest && $(MAKE) tags
#	-@make tags $(MFLGS) $(SRCDIR)/nicelib
#	-@make tags $(MFLGS) $(SRCDIR)/niced
#	-@make tags $(MFLGS) $(SRCDIR)/ntest

######################################################################
# Distribution tools. Assumes your master source is in a directory
# called "nice-$(NICEVERSION)".
TAR	=  $(foreach ENTRY,Makefile ReadMe config.h .maj .min .pat .dat nicelib etc niced ntest,nice-$(NICEVERSION)/$(ENTRY))

.PHONY: tar
tar: clean
	-@tar -cf nice-$(NICEVERSION).tar -C .. $(TAR)
	-@ls -lsag *-$(NICEVERSION).tar

.PHONY: tgz
tgz: clean
	-@tar -czf nice-$(NICEVERSION).tgz -C .. $(TAR)
	-@ls -lsag *-$(NICEVERSION).tgz

######################################################################
# Windows installation. Must ensure first that systype is cygwin and
# that the system has been made and installed into the appropriate
# PREFIX directory, usually "C:\Program Files\Nice $(NICEVERSION)".
.PHONY: zip
zip: 
	-@zip -j nice-$(NICEVERSION).zip /bin/cygwin1.dll
	-@zip -r nice-$(NICEVERSION).zip $(PREFIX)/bin $(PREFIX)/sbin $(PREFIX)/nice $(PREFIX)/etc

######################################################################
# Removes installation.
.PHONY: remove
remove: clean
	-@rm -f $(INCDIR)/niceaux.h
	-@rm -f $(INCDIR)/nicecom.h
	-@rm -f $(INCDIR)/niceapi.h
	-@if [ -d $(INCDIR) ] ; then rmdir $(RDFLGS) $(INCDIR); fi
	-cd $(SRCDIR)/nicelib && $(MAKE) remove
	#-@make $(MFLGS) $(SRCDIR)/nicelib remove
	-@if [ -d $(LIBDIR0) ] ; then rmdir $(RDFLGS) $(LIBDIR0); fi
	-@rm -f $(LIBDIR)
	-cd $(SRCDIR)/niced && $(MAKE) remove
	#-@make $(MFLGS) $(SRCDIR)/niced remove
	-cd $(SRCDIR)/ntest && $(MAKE) remove
	#-@make $(MFLGS) $(SRCDIR)/ntest remove
	-@if [ -d $(BINDIR) ] ; then rmdir $(RDFLGS) $(BINDIR); fi
	-@if [ -d $(SBINDIR) ] ; then rmdir $(RDFLGS) $(SBINDIR); fi
	-@if [ -d $(BINDIR) ] ; then rmdir $(RDFLGS) $(BINDIR); fi
	-@if [ -d $(SBINDIR) ] ; then rmdir $(RDFLGS) $(SBINDIR); fi
	-@if [ -d $(NBINDIR) ] ; then rmdir $(RDFLGS) $(NBINDIR); fi
	-@if [ -d $(NCFGDIR) ] ; then rmdir $(RDFLGS) $(NCFGDIR); fi
	-@if [ -d $(PREFIX) ] ; then rmdir $(RDFLGS) $(PREFIX); fi
ifeq ($(SYSTYPE),debian)
	-@if [ -f /etc/init.d/niced ] ; then rm /etc/init.d/niced ; fi
	-@/usr/sbin/update-rc.d niced remove
endif
ifeq ($(SYSTYPE), redhat)
	-@if [ -f /etc/init.d/niced ] ; then rm /etc/init.d/niced ; fi
	-@rm /etc/rc?.d/???niced 
endif
	-@if [ -f $(CFGDIR)/niced.conf ] ; then echo -e "\nNot removing $(CFGDIR)/niced.conf\n"; fi
	-@if [ -d $(CFGDIR) ] ; then rmdir $(RDFLGS) $(CFGDIR); fi
	-@if [ -d $(NICEDIR) ] ; then rmdir $(RDFLGS) $(NICEDIR); fi
ifeq ($(SYSTYPE), debug)
else
	-@if grep $(NICEUSR) /etc/passwd > /dev/null ; then userdel $(NICEUSR) ; fi
	-@if grep $(NICEGRP) /etc/group > /dev/null ; then groupdel $(NICEGRP) ; fi
endif

######################################################################
# Produces 2up postscript of source.
.PHONY: ps
ps: clean
	-@enscript -Ec -2r -o nice-`date +"%Y%m%d"`.ps $(SRCDIR)/config.h $(SRCDIR)/nicelib/{Makefile,*.[ch]} $(SRCDIR)/niced/*.c $(SRCDIR)/ntest/{Makefile,ntest.c} 

# Produces HTML copy of source
.PHONY: html
html: clean
	-@enscript -Ec -Whtml -o nice-`date +"%Y%m%d"`.html $(SRCDIR)/config.h $(SRCDIR)/nicelib/{Makefile,*.[ch]} $(SRCDIR)/niced/*.c $(SRCDIR)/ntest/{Makefile,ntest.c} 

######################################################################
# Creates a new release.
.PHONY: major
major:
	-@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NICE).1
	-@cat .maj > /tmp/$(NICE).2
	-@echo "1 + p" >> /tmp/$(NICE).2
	-@dc < /tmp/$(NICE).2 > .maj
	-@echo "0" > .min
	-@echo "0" > .pat
	-@date > .dat
	-@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NICE).2
	-@echo "Version" `cat /tmp/$(NICE).2`";" `cat .dat`
	-@mv ../$(NICE)-`cat /tmp/$(NICE).1` ../$(NICE)-`cat /tmp/$(NICE).2`
	-@cd ../$(NICE)-`cat /tmp/$(NICE).2`
	-@sed -e"s%`cat /tmp/$(NICE).1`%`cat /tmp/$(NICE).2`%g" < Makefile > /tmp/$(NICE).3
	-@mv /tmp/$(NICE).3 Makefile
	-@rm /tmp/$(NICE).1 /tmp/$(NICE).2

# Creates a new minor version.
.PHONY: minor
minor:
	-@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NICE).1
	-@cat .min > /tmp/$(NICE).2
	-@echo "1 + p" >> /tmp/$(NICE).2
	-@dc < /tmp/$(NICE).2 > .min
	-@echo "0" > .pat
	-@date > .dat
	-@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NICE).2
	-@echo "Version" `cat /tmp/$(NICE).2`";" `cat .dat`
	-@mv ../$(NICE)-`cat /tmp/$(NICE).1` ../$(NICE)-`cat /tmp/$(NICE).2`
	-@cd ../$(NICE)-`cat /tmp/$(NICE).2`
	-@sed -e"s%`cat /tmp/$(NICE).1`%`cat /tmp/$(NICE).2`%g" < Makefile > /tmp/$(NICE).3
	-@mv /tmp/$(NICE).3 Makefile
	-@rm /tmp/$(NICE).1 /tmp/$(NICE).2

# Creates a new patch version.
.PHONY: patch
patch:
	-@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NICE).1
	-@cat .pat > /tmp/$(NICE).2
	-@echo "1 + p" >> /tmp/$(NICE).2
	-@dc < /tmp/$(NICE).2 > .pat
	-@date > .dat
	-@echo `cat .maj`.`cat .min`.`cat .pat` > /tmp/$(NICE).2
	-@echo "Version" `cat /tmp/$(NICE).2`";" `cat .dat`
	-@mv ../$(NICE)-`cat /tmp/$(NICE).1` ../$(NICE)-`cat /tmp/$(NICE).2`
	-@cd ../$(NICE)-`cat /tmp/$(NICE).2`
	-@sed -e"s%`cat /tmp/$(NICE).1`%`cat /tmp/$(NICE).2`%g" < Makefile > /tmp/$(NICE).3
	-@mv /tmp/$(NICE).3 Makefile
	-@rm /tmp/$(NICE).1 /tmp/$(NICE).2

# Check variable values.
.PHONY: check
check:
	@echo "Version $(NICEVERSION)-$(NICEARCH) for $(NICEOS)"
	@echo "$(NICEDATE)"
	@echo MFLGS=$(MFLGS)
	@echo SYSTYPE=$(SYSTYPE)
	@echo DISTDIR=$(DISTDIR)
	@echo SRCDIR=$(SRCDIR)
	@echo NICEDIR=$(NICEDIR)
	@echo SBINDIR=$(SBINDIR)
	@echo BINDIR=$(BINDIR)
	@echo CFGDIR=$(CFGDIR)
	@echo INCDIR=$(INCDIR)
	@echo LIBDIR=$(LIBDIR)
	@echo NBINDIR=$(NBINDIR)
	@echo NCFGDIR=$(NCFGDIR)
