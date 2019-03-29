# Kelvin Utilities
# Alberto Maria Segre
#
# Copyright (c) 2005, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

######################################################################
# No user changes below this point!

# Object files
OBJ = utils.o

# Include files
INC = utils.h

# Libraries
LIB =

######################################################################
# Make an object file from source
%.o      : %.c $(INCS)
	-@$(CC) -c $(CFLGS) -o $@ $<

######################################################################
.PHONY: install
install: libutils.a
	-@cp utils.h $(INCDIR)
	-@cp libutils.a $(LIBDIR)

libutils.a: Makefile $(OBJ)
	-@ar -rsc libutils.a $(OBJ)

# To ensure a clean build, use "make fresh."
.PHONY: fresh
fresh: clean install

# Deletes all compiled files.
.PHONY: clean
clean:
	-@rm -f libutils.a *.o *~ 
	-@rm -f $(INCDIR)/utils.h
	-@rm -f $(LIBDIR)/libutils.a