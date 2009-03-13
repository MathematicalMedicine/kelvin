# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = config.h model_options.h model_type.h model_range.h lambda_cell.h

OBJS = config.o read_pedfile.o

all : $(KVNLIBDIR)/libconfig.a $(KVNINCDIR)/config.h \
      $(KVNINCDIR)/lambda_cell.h $(KVNINCDIR)/model_range.h $(KVNINCDIR)/model_type.h \
      $(KVNINCDIR)/model_options.h

clean :
	rm -f $(OBJS) libconfig.a $(KVNLIBDIR)/libconfig.a


$(KVNINCDIR)/config.h : config.h
	 cp config.h $(KVNINCDIR)

$(KVNINCDIR)/lambda_cell.h : lambda_cell.h
	 cp lambda_cell.h $(KVNINCDIR)

$(KVNINCDIR)/model_range.h : model_range.h
	 cp model_range.h $(KVNINCDIR)

$(KVNINCDIR)/model_type.h : model_type.h
	 cp model_type.h $(KVNINCDIR)

$(KVNINCDIR)/model_options.h : model_options.h
	 cp model_options.h $(KVNINCDIR)

$(KVNLIBDIR)/libconfig.a : libconfig.a
	cp libconfig.a $(KVNLIBDIR)

libconfig.a : $(OBJS)
	ar -rsc libconfig.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(EXTRAFLAG) $< -o $@
