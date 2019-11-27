# NICE System
# Alberto Maria Segre
#
# Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

BIN	= ntest
LIB	= -lniceapi -lnicecom -lniceaux
SRC	= ntest.c
INC	= 
TAR	= Makefile $(SRC) $(INC)

######################################################################
# Determine application name. This is the compiled code's notion of
# who I am. Normally, $(NICEVERSION) and $(NICEARCH) would come from
# the application Makefile (some sort of application versioning system
# required, and uname -m can be used to retrieve the
# architecture). But for ntest, which is part of the distribution,
# we'll just use inherited versions from the infrastructure, meaning
# ntest versioning is locked in step with the infraastructure version.

NICEAPP = $(BIN)-$(NICEVERSION)-$(NICEARCH)
ifneq ($(NICEOS),Cygwin)
  NICEEXE = $(NBINDIR)/$(NICEAPP)
else
  NICEEXE = $(NBINDIR)/$(NICEAPP).exe
endif

# NOTE: I (Jo) have disabled this default target, as "install" should not be
# the default.
#######################################################################
## Default target.  Installs application in $NBINDIR.
#install: ntest
#	@cp $(NICEAPP) $(NICEEXE)
#	@chmod 755 $(NICEEXE)

######################################################################
# Sample NICE application. Note that $(NICEAPP) is not passed into the
# application itself. Everything is really driven off of the executable
# filename on the root processor.
ntest: $(INC) $(SRC)
	@$(CC) $(CFLGS) $(SRC) -o $(NICEAPP) $(LIB)

######################################################################
# Remove installed version; leaves local binary.
remove:
	-@rm $(NICEEXE)

######################################################################
# Utilities.
clean:
	-@rm -f $(NICEAPP) TAGS *~ ntest-????????.ps

ps:
	@enscript -Ec -2r -o ntest-`date +"%Y%m%d"`.ps $(TAR)

tags:
	@etags *.[ch]
