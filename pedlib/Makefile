# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = pedigree.h locus.h genotype_elimination.h \
       likelihood.h pedlib.h allele_set_recoding.h

OBJS = locus.o allele_set_recoding.o \
       genotype_elimination.o parental_pair.o likelihood.o

all : $(KVNLIBDIR)/libped.a $(KVNINCDIR)/pedigree.h $(KVNINCDIR)/locus.h \
      $(KVNINCDIR)/genotype_elimination.h \
      $(KVNINCDIR)/likelihood.h $(KVNINCDIR)/pedlib.h \
      $(KVNINCDIR)/allele_set_recoding.h

clean :
	rm -f $(OBJS) libped.a $(KVNLIBDIR)/libped.a

$(KVNINCDIR)/pedigree.h : pedigree.h
	 cp pedigree.h $(KVNINCDIR)

$(KVNINCDIR)/locus.h : locus.h
	 cp locus.h $(KVNINCDIR)

$(KVNINCDIR)/genotype_elimination.h : genotype_elimination.h
	 cp genotype_elimination.h $(KVNINCDIR)

$(KVNINCDIR)/likelihood.h : likelihood.h
	 cp likelihood.h $(KVNINCDIR)

$(KVNINCDIR)/pedlib.h : pedlib.h
	 cp pedlib.h $(KVNINCDIR)

$(KVNINCDIR)/allele_set_recoding.h : allele_set_recoding.h
	 cp allele_set_recoding.h $(KVNINCDIR)

$(KVNLIBDIR)/libped.a : libped.a
	cp libped.a $(KVNLIBDIR)

libped.a : $(OBJS)
	ar -rsc libped.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) $< -o $@
