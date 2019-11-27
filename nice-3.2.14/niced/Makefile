# NICE System
# Alberto Maria Segre
#
# Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

LIB	= -lnicecom -lniceaux
APILIB	= -lniceapi
SRC	= niced.c nasti.c bench.c niceq.c
INC	= niced.h
ifneq ($(OS),Cygwin)
SBIN	= niced
BIN	= niceq
else
SBIN	= niced.exe
BIN	= niceq.exe
endif
TAR	= Makefile $(SRC) $(INC)
#CFLGS   += -DNOFORK
# Location of NICE-enabled applications.
CFLGS	+= -DNBINDIR=\"$(NBINDIR)\" -DCFGDIR=\"$(CFGDIR)\"
# Version and architecture information.
CFLGS	+= -DNICEVERSION=\"$(NICEVERSION)\" -DNICEARCH=\"$(NICEARCH)\"
# NOTE: I (Jo) have disabled this default target, as "install" should not be
# the default.
######################################################################
# Default target.  Installs NICE daemon in $SBINDIR, typically
# /usr/local/sbin.  Installs NICE query program in $BINDIR,
# typically /usr/local/bin. 
#.PHONY: all
#all: $(SBIN) $(BIN)
#
#.PHONY: install
#install: $(SBIN) $(BIN)
#	@cp $(SBIN) $(SBINDIR)
#	@(cd $(SBINDIR))
#	@cp $(BIN) $(BINDIR)
#	@(cd $(BINDIR))
#
######################################################################
# The NICE daemon and the NICE query program.
niced niced.exe niceq niceq.exe: niced.h niced.c niceq.c nasti.c bench.c
	@$(CC) $(CFLGS) niced.c nasti.c bench.c -o niced $(LIB)
	@$(CC) $(CFLGS) niceq.c -o niceq $(LIB)

######################################################################
# Remove installed daemon and query program.
.PHONY: remove
remove:
	-@(cd $(SBINDIR); rm -f $(SBIN))
	-@(cd $(BINDIR); rm -f $(BIN))

######################################################################
# Utilities.
.PHONY: clean
clean:
	-@rm -f $(SBIN) $(BIN) TAGS *~ niced-????????.ps

.PHONY: ps
ps:
	@enscript -Ec -2r -o niced-`date +"%Y%m%d"`.ps $(TAR)

.PHONY: tags
tags:
	@etags *.[ch]
