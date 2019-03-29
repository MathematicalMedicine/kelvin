# Pedigree Library
# Yungui Huang
#
# Copyright (c) 2005, The University of Iowa.  All rights reserved.
# Permission is hereby given to use and reproduce this software 
# for non-profit educational purposes only.

# Compile options:
#
# -DNO_POLYNOMIAL removes the polynomial version of the code.
##CC = gcc
#CFLGS = -g -Wall
#CFLGS += -DNO_POLYNOMIAL

######################################################################
# No user changes below this point!

# Object files
OBJ = read_pedfile.o locus.o tools.o allele_set_recoding.o \
	genotype_elimination.o parental_pair.o likelihood.o polynomial.o

# Include files
INCS = pedigree.h locus.h model_options.h  \
	genotype_elimination.h likelihood.h polynomial.h \
	pedlib.h allele_set_recoding.h tools.h 

# Libraries
LIB = -lutils -lped -lgsl -lgslcblas -lm 

######################################################################
# Make an object file from source
%.o : %.c $(INCS)
	-@$(CC) -c $(CFLGS) -I ../utils -o $@ $<

######################################################################
# Default target
.PHONY: install
install: libped.a
	-@cp {pedlib,pedigree,locus,model_options,polynomial,likelihood}.h $(INCDIR)
	-@cp libped.a $(LIBDIR)

libped.a: Makefile $(OBJ)
	-@ar -rsc libped.a $(OBJ)


######################################################################
# Utilities.

# To ensure a clean build, use "make fresh."
.PHONY: fresh
fresh: clean libped.a

# Deletes all compiled files.
.PHONY: clean
clean:
	-@rm -f libped.a *.o *~ test
	-@rm -f $(INCDIR)/{pedlib,pedigree,locus,model_options,likelihood,polynomial}.h
	-@rm -f $(LIBDIR)/libped.a
