# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

INCS = pedigree.h locus.h genotype_elimination.h \
       likelihood.h pedlib.h allele_set_recoding.h

OBJS = locus.o allele_set_recoding.o \
       genotype_elimination.o parental_pair.o likelihood.o

all : ../lib/libped.a

clean :
	rm -f $(OBJS) libped.a ../lib/libped.a

../lib/libped.a : libped.a
	cp libped.a ../lib

libped.a : $(OBJS)
	ar -rsc libped.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) $< -o $@
