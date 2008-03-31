# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = pedigree.h locus.h model_options.h genotype_elimination.h \
       likelihood.h polynomial.h pedlib.h allele_set_recoding.h tools.h

OBJS = read_pedfile.o locus.o tools.o allele_set_recoding.o \
       genotype_elimination.o parental_pair.o likelihood.o polynomial.o

all : $(KVNLIBDIR)/libped.a $(KVNINCDIR)/pedigree.h $(KVNINCDIR)/locus.h \
      $(KVNINCDIR)/model_options.h $(KVNINCDIR)/genotype_elimination.h \
      $(KVNINCDIR)/likelihood.h $(KVNINCDIR)/polynomial.h $(KVNINCDIR)/pedlib.h \
      $(KVNINCDIR)/allele_set_recoding.h $(KVNINCDIR)/tools.h

clean :
	rm -f $(OBJS) libped.a $(KVNLIBDIR)/libped.a

$(KVNINCDIR)/pedigree.h : pedigree.h
	 cp pedigree.h $(KVNINCDIR)

$(KVNINCDIR)/locus.h : locus.h
	 cp locus.h $(KVNINCDIR)

$(KVNINCDIR)/model_options.h : model_options.h
	 cp model_options.h $(KVNINCDIR)

$(KVNINCDIR)/genotype_elimination.h : genotype_elimination.h
	 cp genotype_elimination.h $(KVNINCDIR)

$(KVNINCDIR)/likelihood.h : likelihood.h
	 cp likelihood.h $(KVNINCDIR)

$(KVNINCDIR)/polynomial.h : polynomial.h
	 cp polynomial.h $(KVNINCDIR)

$(KVNINCDIR)/pedlib.h : pedlib.h
	 cp pedlib.h $(KVNINCDIR)

$(KVNINCDIR)/allele_set_recoding.h : allele_set_recoding.h
	 cp allele_set_recoding.h $(KVNINCDIR)

$(KVNINCDIR)/tools.h : tools.h
	 cp tools.h $(KVNINCDIR)

$(KVNLIBDIR)/libped.a : libped.a
	cp libped.a $(KVNLIBDIR)

libped.a : $(OBJS)
	ar -rsc libped.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@
