# NICE System
# Alberto Maria Segre
#
# Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

SRC	= niceaux.c nicecom.c niceapi.c
INC	= niceaux.h nicecom.h niceapi.h 
TAR	= Makefile $(SRC) $(INC)

# __USE_ISOC99  may be needed for vsscanf definition in stdio.h?
#
CFLGS	+= -D__USE_ISOC99

######################################################################
# NOTE: I (Jo) have disabled this default target, as "install" should not be
# the default.
## Default target.
#.PHONY: install
#install: nicelib
#	@cp libnicecom.a libniceaux.a libniceapi.a $(LIBDIR)

.PHONY: nicelib
nicelib: libnicecom.a libniceaux.a libniceapi.a

.PHONY: remove
remove:
	-@(cd $(LIBDIR); rm -f libnicecom.a libniceaux.a libniceapi.a)

######################################################################
# Libraries.
libnicecom.a: nicecom.o
	@ar -rsc libnicecom.a nicecom.o

libniceaux.a: niceaux.o
	@ar -rsc libniceaux.a niceaux.o

libniceapi.a: niceapi.o
	@ar -rsc libniceapi.a niceapi.o 

# NICE communication library.
nicecom.o: nicecom.c 
	@$(CC) -c $(CFLGS) nicecom.c

# NICE software aux library.
niceaux.o: niceaux.c 
	@$(CC) -c $(CFLGS) niceaux.c

# NICE API.
niceapi.o: nicecom.c niceapi.c 
	@$(CC) -c $(CFLGS) niceapi.c 

######################################################################
# Utilities.
.PHONY: clean
clean:
	-@rm -f *~ nicecom.o niceaux.o niceapi.o TAGS
	-@rm -f libnicecom.a libniceaux.a libniceapi.a nicelib-????????.ps

.PHONY: ps
ps:
	@enscript -Ec -2r -o nicelib-`date +"%Y%m%d"`.ps Makefile $(INC) $(SRC)
#	@enscript -Ec -o nicelib-`date +"%Y%m%d"`.ps Makefile $(INC) $(SRC)

.PHONY: tags
tags:
	@etags *.[ch]
