# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = utils.h sw.h tpl.h hashtab.h lookupa.h recycle.h standard.h sSDHandler.h dists.h
OBJS = utils.o sw.o tpl.o hashtab.o lookupa.o recycle.o sSDHandler.o dcdflib.o dists.o ipmpar.o

all : $(KVNLIBDIR)/libutils.a \
      $(KVNINCDIR)/utils.h \
      $(KVNINCDIR)/sw.h \
      $(KVNINCDIR)/sSDHandler.h \
      $(KVNINCDIR)/dists.h \
      $(KVNINCDIR)/tpl.h \
      $(KVNINCDIR)/hashtab.h \
      $(KVNINCDIR)/lookupa.h \
      $(KVNINCDIR)/recycle.h \
      $(KVNINCDIR)/standard.h

$(KVNLIBDIR)/libutils.a : libutils.a
	cp libutils.a $(KVNLIBDIR)

$(KVNINCDIR)/utils.h : utils.h
	cp utils.h $(KVNINCDIR)

$(KVNINCDIR)/sw.h : sw.h
	cp sw.h $(KVNINCDIR)

$(KVNINCDIR)/sSDHandler.h : sSDHandler.h
	cp sSDHandler.h $(KVNINCDIR)

$(KVNINCDIR)/dists.h : dists.h
	cp dists.h $(KVNINCDIR)

$(KVNINCDIR)/tpl.h : tpl.h
	cp tpl.h $(KVNINCDIR)

$(KVNINCDIR)/hashtab.h : hashtab.h
	cp hashtab.h $(KVNINCDIR)

$(KVNINCDIR)/lookupa.h : lookupa.h
	cp lookupa.h $(KVNINCDIR)

$(KVNINCDIR)/recycle.h : recycle.h
	cp recycle.h $(KVNINCDIR)

$(KVNINCDIR)/standard.h : standard.h
	cp standard.h $(KVNINCDIR)

libutils.a : $(OBJS)
	ar -rsc libutils.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -f $(OBJS) libutils.a $(KVNLIBDIR)/libutils.a
