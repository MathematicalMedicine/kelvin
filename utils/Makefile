# Utilities
# Alberto Maria Segre et. al.
#
# Copyright (c) 2005, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

######################################################################
# No user changes below this point!

# Object files
OBJ = utils.o sw.o tpl.o hashtab.o lookupa.o recycle.o

# Include files
INC = utils.h sw.h tpl.h hashtab.h lookupa.h recycle.h standard.h

# Libraries
LIB =

######################################################################
# Make an object file from source
%.o      : %.c $(INC)
	-@$(CC) -c $(CFLGS) -o $@ $<

######################################################################
.PHONY: install
install: libutils.a
	-@cp $(INC) $(INCDIR)
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
	-@rm -f $(LIBDIR)/libutils.a